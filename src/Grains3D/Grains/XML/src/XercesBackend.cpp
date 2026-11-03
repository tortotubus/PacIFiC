#include "XercesBackend.hh"
#include "XercesRuntime.hh"

#include <xercesc/dom/DOM.hpp>
#include <xercesc/framework/MemBufInputSource.hpp>
#include <xercesc/framework/MemBufFormatTarget.hpp>
#include <xercesc/internal/XMLScanner.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/sax/ErrorHandler.hpp>
#include <xercesc/sax/Locator.hpp>
#include <xercesc/sax/SAXParseException.hpp>
#include <xercesc/util/TransService.hpp>
#include <xercesc/util/XMLChar.hpp>
#include <xercesc/util/XMLResourceIdentifier.hpp>
#include <xercesc/util/XMLEntityResolver.hpp>
#include <xercesc/util/XMLUri.hpp>
#include <xercesc/util/SecurityManager.hpp>
#include <xercesc/validators/DTD/DTDEntityDecl.hpp>

#include <fstream>
#include <iterator>
#include <unordered_map>

namespace grains::xml {
namespace {
using namespace xercesc;
using Runtime = grains_xml_internal::XercesRuntime;

template<class T> struct Release {
    void operator()(T* value) const noexcept { if (value) value->release(); }
};
template<class T> using Owned = std::unique_ptr<T, Release<T>>;

std::string utf8(const XMLCh* value) {
    if (!value) return {};
    TranscodeToStr converted(value, "UTF-8");
    return {reinterpret_cast<const char*>(converted.str()), converted.length()};
}

std::shared_ptr<Runtime> makeRuntime() {
    try { return std::make_shared<Runtime>(); }
    catch (const XMLException&) { throw Error("Cannot initialize Xerces runtime"); }
}

// Explicit UTF-8 conversion, independent of the host locale. Validate before
// passing strings to DOM APIs, which do not enforce all XML character rules.
std::basic_string<XMLCh> wide(std::string_view value) {
    if (value.find('\0') != std::string_view::npos) throw Error("Invalid XML character");
    try {
        TranscodeFromStr converted(reinterpret_cast<const XMLByte*>(value.data()), value.size(), "UTF-8");
        std::basic_string<XMLCh> result(converted.str(), converted.length());
        for (std::size_t i = 0; i < result.size(); ++i) {
            const auto ch = static_cast<unsigned int>(result[i]);
            if (ch >= 0xD800 && ch <= 0xDBFF && i + 1 < result.size() &&
                result[i + 1] >= 0xDC00 && result[i + 1] <= 0xDFFF) { ++i; continue; }
            if (!(ch == 9 || ch == 10 || ch == 13 || (ch >= 0x20 && ch <= 0xD7FF) ||
                  (ch >= 0xE000 && ch <= 0xFFFD))) throw Error("Invalid XML character");
        }
        return result;
    } catch (const XMLException&) { throw Error("Invalid UTF-8 XML text"); }
}
std::basic_string<XMLCh> name(std::string_view value) {
    auto result = wide(value);
    if (result.empty() || value.find(':') != std::string_view::npos || value == "xmlns" ||
        !XMLChar1_0::isValidNCName(result.c_str(), result.size()))
        throw Error("Invalid or unsupported XML name: " + std::string(value));
    return result;
}

std::string readFile(const std::filesystem::path& path) {
    std::error_code error;
    if (!std::filesystem::is_regular_file(path, error))
        throw Error("Cannot open XML resource as a regular file", {path.string(), {}, {}});
    std::ifstream input(path, std::ios::binary);
    if (!input) throw Error("Cannot open XML resource", {path.string(), {}, {}});
    std::string data{std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};
    if (input.bad()) throw Error("Cannot read XML resource", {path.string(), {}, {}});
    return data;
}
std::string fileUri(const std::filesystem::path& path) {
    const char hex[] = "0123456789ABCDEF";
    std::string result = "file://";
    for (unsigned char c : path.generic_u8string()) {
        if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
            (c >= '0' && c <= '9') || c == '/' || c == ':' || c == '-' || c == '_' || c == '.' || c == '~')
            result += static_cast<char>(c);
        else { result += '%'; result += hex[c >> 4]; result += hex[c & 15]; }
    }
    return result;
}
std::filesystem::path resourcePath(const XMLCh* system, const XMLCh* base) {
    if (!system) throw Error("Missing XML resource identifier");
    XMLUri baseUri(base);
    XMLUri uri(&baseUri, system);
    if (utf8(uri.getScheme()) != "file" ||
        (!utf8(uri.getHost()).empty() && utf8(uri.getHost()) != "localhost") ||
        uri.getUserInfo() || uri.getQueryString() || uri.getFragment() || uri.getPort() != -1)
        throw Error("XML resource is not permitted");
    auto path = utf8(uri.getPath());
    std::string decoded;
    const auto digit = [](char c) -> int {
        if (c >= '0' && c <= '9') return c - '0';
        if (c >= 'a' && c <= 'f') return c - 'a' + 10;
        if (c >= 'A' && c <= 'F') return c - 'A' + 10;
        return -1;
    };
    for (std::size_t i = 0; i < path.size(); ++i) {
        if (path[i] != '%') decoded += path[i];
        else {
            if (i + 2 >= path.size() || digit(path[i+1]) < 0 || digit(path[i+2]) < 0)
                throw Error("Invalid resource URI");
            const char c = static_cast<char>(digit(path[i+1]) * 16 + digit(path[i+2]));
            if (!c) throw Error("Invalid resource URI");
            decoded += c;
            i += 2;
        }
    }
    return std::filesystem::u8path(decoded);
}

class Parser final : public XercesDOMParser, public ErrorHandler, public XMLEntityResolver {
public:
    std::filesystem::path source;
    std::optional<std::filesystem::path> dtdRoot;
    std::unordered_map<const DOMNode*, std::size_t> lines;
    Parser() { setErrorHandler(this); setXMLEntityResolver(this); setDisableDefaultEntityResolution(true); }
    void warning(const SAXParseException&) override {}
    void error(const SAXParseException& e) override { fail(e); }
    void fatalError(const SAXParseException& e) override { fail(e); }
    void resetErrors() override {}
    [[noreturn]] void fail(const SAXParseException& e) {
        throw Error(utf8(e.getMessage()), {source.string(), {}, static_cast<std::size_t>(e.getLineNumber())});
    }
    void entityDecl(const DTDEntityDecl& declaration, bool parameter, bool ignored) override {
        // Reject external general declarations before the scanner can resolve
        // them; Xerces's resource callback does not distinguish general/PEs.
        if (!parameter && !ignored && declaration.getSystemId())
            throw Error("External general XML entities are not permitted", {source.string(), {}, {}});
        XercesDOMParser::entityDecl(declaration, parameter, ignored);
    }
    InputSource* resolveEntity(XMLResourceIdentifier* resource) override {
        if (!dtdRoot || resource->getResourceIdentifierType() != XMLResourceIdentifier::ExternalEntity)
            throw Error("XML resource is not permitted", {source.string(), {}, {}});
        std::filesystem::path path;
        try {
            const auto fallback = wide(fileUri(source));
            path = std::filesystem::canonical(resourcePath(resource->getSystemId(),
                resource->getBaseURI() ? resource->getBaseURI() : fallback.c_str()));
        } catch (const Error&) { throw; }
        catch (...) { throw Error("Cannot open XML resource", {source.string(), {}, {}}); }
        const auto relative = path.lexically_relative(*dtdRoot);
        if (relative.empty() || relative.is_absolute() || *relative.begin() == "..")
            throw Error("XML resource is not permitted", {path.string(), {}, {}});
        const auto data = readFile(path);
        const auto uri = wide(fileUri(path));
        // The input source adopts its own buffer, independent of this callback.
        auto buffer = std::make_unique<XMLByte[]>(data.size());
        std::copy(data.begin(), data.end(), buffer.get());
        auto input = std::make_unique<MemBufInputSource>(buffer.get(), data.size(), uri.c_str(), true);
        buffer.release();
        return input.release();
    }
    void startElement(const XMLElementDecl& decl, unsigned int url,
                      const XMLCh* const prefix, const RefVectorOf<XMLAttr>& attrs,
                      XMLSize_t count, bool empty, bool root) override {
        XercesDOMParser::startElement(decl, url, prefix, attrs, count, empty, root);
        const auto* node = getCurrentNode();
        lines[node] = static_cast<std::size_t>(getScanner()->getLocator()->getLineNumber());
    }
};

class Storage final : public detail::WritableStorage {
public:
    Storage(std::shared_ptr<Runtime> runtime, Owned<DOMDocument> document,
            std::filesystem::path source = {}, std::unordered_map<const DOMNode*, std::size_t> lines = {})
        : runtime_(std::move(runtime)), document_(std::move(document)), source_(std::move(source)), lines_(std::move(lines)) {
        std::vector<DOMNode*> pending{document_->getDocumentElement()};
        while (!pending.empty()) {
            auto* node = pending.back(); pending.pop_back();
            ids_.emplace(node, nodes_.size()); nodes_.push_back(node);
            for (auto* child = node->getLastChild(); child; child = child->getPreviousSibling())
                if (child->getNodeType() == DOMNode::ELEMENT_NODE) pending.push_back(child);
        }
    }
    detail::NodeId root() const override { return 0; }
    std::string name(detail::NodeId id) const override { return utf8(nodes_.at(id)->getNodeName()); }
    std::optional<std::string> attribute(detail::NodeId id, std::string_view key) const override {
        const auto attr = wide(key);
        const auto* node = nodes_.at(id)->getAttributes()->getNamedItem(attr.c_str());
        if (!node) return std::nullopt;
        return utf8(node->getNodeValue());
    }
    std::string text(detail::NodeId id) const override {
        std::string result;
        for (auto* child = nodes_.at(id)->getFirstChild(); child; child = child->getNextSibling())
            if (child->getNodeType() == DOMNode::TEXT_NODE || child->getNodeType() == DOMNode::CDATA_SECTION_NODE)
                result += utf8(child->getNodeValue());
        return result;
    }
    std::vector<detail::NodeId> children(detail::NodeId id) const override {
        std::vector<detail::NodeId> result;
        for (auto* child = nodes_.at(id)->getFirstChild(); child; child = child->getNextSibling())
            if (child->getNodeType() == DOMNode::ELEMENT_NODE) result.push_back(ids_.at(child));
        return result;
    }
    std::optional<detail::NodeId> parent(detail::NodeId id) const override {
        const auto found = ids_.find(nodes_.at(id)->getParentNode());
        return found == ids_.end() ? std::nullopt : std::optional<detail::NodeId>(found->second);
    }
    std::optional<detail::NodeId> firstChild(detail::NodeId id) const override {
        for (auto* node = nodes_.at(id)->getFirstChild(); node; node = node->getNextSibling())
            if (node->getNodeType() == DOMNode::ELEMENT_NODE) return ids_.at(node);
        return std::nullopt;
    }
    std::optional<detail::NodeId> nextSibling(detail::NodeId id) const override {
        for (auto* node = nodes_.at(id)->getNextSibling(); node; node = node->getNextSibling())
            if (node->getNodeType() == DOMNode::ELEMENT_NODE) return ids_.at(node);
        return std::nullopt;
    }
    SourceLocation location(detail::NodeId id) const override {
        const auto* node = nodes_.at(id);
        std::string path;
        for (const auto* current = node; current && current->getNodeType() == DOMNode::ELEMENT_NODE;
             current = current->getParentNode()) {
            std::size_t index = 1;
            for (auto* sibling = current->getPreviousSibling(); sibling; sibling = sibling->getPreviousSibling())
                if (sibling->getNodeType() == DOMNode::ELEMENT_NODE &&
                    XMLString::equals(sibling->getNodeName(), current->getNodeName())) ++index;
            path = "/" + utf8(current->getNodeName()) + "[" + std::to_string(index) + "]" + path;
        }
        SourceLocation result{source_.string(), path, {}};
        const auto found = lines_.find(node);
        if (found != lines_.end()) result.line = found->second;
        return result;
    }
    detail::NodeId appendChild(detail::NodeId id, std::string_view key) override {
        const auto valid = ::grains::xml::name(key);
        auto* parent = nodes_.at(id);
        Owned<DOMElement> node(document_->createElement(valid.c_str()));
        const auto newId = nodes_.size();
        nodes_.push_back(node.get());
        try { ids_.emplace(node.get(), newId); }
        catch (...) { nodes_.pop_back(); throw; }
        try { parent->appendChild(node.get()); }
        catch (...) { ids_.erase(node.get()); nodes_.pop_back(); throw; }
        node.release();
        return newId;
    }
    void setAttribute(detail::NodeId id, std::string_view key, std::string_view value) override {
        const auto valid = ::grains::xml::name(key), text = wide(value);
        static_cast<DOMElement*>(nodes_.at(id))->setAttribute(valid.c_str(), text.c_str());
    }
    void appendText(detail::NodeId id, std::string_view value) override {
        const auto text = wide(value);
        auto* parent = nodes_.at(id);
        Owned<DOMText> node(document_->createTextNode(text.c_str()));
        parent->appendChild(node.get());
        node.release();
    }
    std::string serialize(bool pretty) const override {
        try {
            const auto ls = wide("LS");
            auto* implementation = dynamic_cast<DOMImplementationLS*>(DOMImplementationRegistry::getDOMImplementation(ls.c_str()));
            if (!implementation) throw Error("Xerces DOM LS is unavailable");
            Owned<DOMLSSerializer> serializer(implementation->createLSSerializer());
            // Xerces pretty-printing may change mixed content. Compact output
            // is permitted by the API even when formatting is requested.
            (void)pretty;
            serializer->getDomConfig()->setParameter(XMLUni::fgDOMWRTFormatPrettyPrint, false);
            Owned<DOMLSOutput> output(implementation->createLSOutput());
            const auto encoding = wide("UTF-8");
            output->setEncoding(encoding.c_str());
            MemBufFormatTarget target;
            output->setByteStream(&target);
            if (!serializer->write(document_.get(), output.get())) throw Error("Cannot serialize XML document");
            return {reinterpret_cast<const char*>(target.getRawBuffer()), target.getLen()};
        } catch (const XMLException& e) { throw Error(utf8(e.getMessage())); }
        catch (const DOMException& e) { throw Error(utf8(e.getMessage())); }
    }
private:
    // Declared first so the native document is released before the runtime.
    std::shared_ptr<Runtime> runtime_;
    Owned<DOMDocument> document_;
    std::filesystem::path source_;
    std::unordered_map<const DOMNode*, std::size_t> lines_;
    std::vector<DOMNode*> nodes_;
    std::unordered_map<const DOMNode*, detail::NodeId> ids_;
};
}

Document XercesBackend::parseFile(const std::filesystem::path& file, const ParseOptions& options) const {
    return parseString(readFile(file), file, options);
}
Document XercesBackend::parseString(std::string_view source, const std::filesystem::path& sourceFile,
                                   const ParseOptions& options) const {
    const auto path = std::filesystem::absolute(sourceFile.empty() ? "memory.xml" : sourceFile).lexically_normal();
    if (options.dtd != DtdPolicy::Disabled && options.dtd != DtdPolicy::Validate && options.dtd != DtdPolicy::LoadDefaults)
        throw Error("Unsupported DTD policy", {path.string(), {}, {}});
    if (options.dtd == DtdPolicy::Disabled && options.dtdDirectory)
        throw Error("dtdDirectory requires DTD loading or validation", {path.string(), {}, {}});
    auto runtime = makeRuntime();
    try {
        SecurityManager security;
        Parser parser;
        parser.source = path;
        if (options.dtdDirectory) {
            std::error_code error;
            auto root = std::filesystem::canonical(*options.dtdDirectory, error);
            if (error || !std::filesystem::is_directory(root)) throw Error("DTD directory does not exist");
            parser.dtdRoot = std::move(root);
        }
        security.setEntityExpansionLimit(10000);
        parser.setSecurityManager(&security);
        parser.setDoNamespaces(true);
        parser.setDoSchema(false);
        parser.setCreateEntityReferenceNodes(false);
        parser.setLoadExternalDTD(options.dtd != DtdPolicy::Disabled);
        parser.setValidationScheme(options.dtd == DtdPolicy::Validate ? XercesDOMParser::Val_Always : XercesDOMParser::Val_Never);
        const auto uri = wide(fileUri(path));
        MemBufInputSource input(reinterpret_cast<const XMLByte*>(source.data()), source.size(), uri.c_str(), false);
        parser.parse(input);
        if (parser.getErrorCount() || !parser.getDocument() || !parser.getDocument()->getDocumentElement())
            throw Error("XML parsing failed", {path.string(), {}, {}});
        return Document(std::make_shared<Storage>(runtime, Owned<DOMDocument>(parser.adoptDocument()), path, std::move(parser.lines)));
    } catch (const XMLException& e) { throw Error(utf8(e.getMessage()), {path.string(), {}, {}}); }
    catch (const DOMException& e) { throw Error(utf8(e.getMessage()), {path.string(), {}, {}}); }
    catch (const SAXException& e) { throw Error(utf8(e.getMessage()), {path.string(), {}, {}}); }
}
WritableDocument XercesBackend::createDocument(std::string_view rootName) const {
    auto runtime = makeRuntime();
    try {
        const auto root = name(rootName), core = wide("Core");
        auto* implementation = DOMImplementationRegistry::getDOMImplementation(core.c_str());
        if (!implementation) throw Error("Xerces DOM Core is unavailable");
        Owned<DOMDocument> document(implementation->createDocument(nullptr, root.c_str(), nullptr));
        return WritableDocument(std::make_shared<Storage>(runtime, std::move(document)));
    } catch (const XMLException& e) { throw Error(utf8(e.getMessage())); }
    catch (const DOMException& e) { throw Error(utf8(e.getMessage())); }
}
} // namespace grains::xml
