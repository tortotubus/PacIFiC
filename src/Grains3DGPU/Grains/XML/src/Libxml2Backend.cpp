#include "Libxml2Backend.hh"

#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/uri.h>
#include <libxml/xmlerror.h>

#include <climits>
#include <exception>
#include <fstream>
#include <iterator>
#include <unordered_map>

#if LIBXML_VERSION < 21400
#error "The Grains XML backend requires libxml2 >= 2.14"
#endif

namespace grains::xml {
namespace {

struct XmlFree {
  void operator()(void *p) const noexcept { xmlFree(p); }
};

using XmlString = std::unique_ptr<xmlChar, XmlFree>;
using NativeDocument = std::unique_ptr<xmlDoc, decltype(&xmlFreeDoc)>;
using NativeContext = std::unique_ptr<xmlParserCtxt, decltype(&xmlFreeParserCtxt)>;
using NativeUri = std::unique_ptr<xmlURI, decltype(&xmlFreeURI)>;

std::string string(const xmlChar *value) {
  return value ? reinterpret_cast<const char *>(value) : "";
}

std::string qualifiedName(const xmlChar *name, const xmlNs *ns) {
  return ns && ns->prefix ? string(ns->prefix) + ":" + string(name)
                          : string(name);
}

std::string readFile(const std::filesystem::path &path) {
  std::error_code error;
  if (!std::filesystem::is_regular_file(path, error))
    throw Error("Cannot open XML resource as a regular file",
                {path.string(), {}, {}});
  std::ifstream input(path, std::ios::binary);
  if (!input)
    throw Error("Cannot open XML resource", {path.string(), {}, {}});
  try {
    std::string contents{std::istreambuf_iterator<char>(input),
                         std::istreambuf_iterator<char>()};
    if (input.bad())
      throw Error("Cannot read XML resource", {path.string(), {}, {}});
    return contents;
  } catch (const std::ios_base::failure &) {
    throw Error("Cannot read XML resource", {path.string(), {}, {}});
  }
}

std::string fileUri(const std::filesystem::path &path) {
  XmlString uri(
      xmlPathToURI(reinterpret_cast<const xmlChar *>(path.u8string().c_str())));
  if (!uri)
    throw Error("Cannot encode resource path", {path.string(), {}, {}});
  return string(uri.get());
}

std::filesystem::path localPath(const char *url) {
  if (!url)
    throw Error("Missing resource URL");
  const auto native = std::filesystem::u8path(url);
  if (native.is_absolute())
    return native;
  NativeUri uri(xmlParseURI(url), xmlFreeURI);
  if (!uri || !uri->path ||
      (uri->scheme && std::string(uri->scheme) != "file") ||
      (uri->server && std::string(uri->server) != "localhost" &&
       *uri->server) ||
      uri->query || uri->fragment)
    throw Error(std::string("Only local file resources are supported: ") + url);
  std::unique_ptr<char, XmlFree> path(
      xmlURIUnescapeString(uri->path, 0, nullptr));
  if (!path)
    throw Error("Cannot decode resource path");
  return std::filesystem::u8path(path.get());
}

struct ParseState {
  std::filesystem::path source;
  std::optional<std::filesystem::path> dtdRoot;
  std::string message;
  SourceLocation location;
  std::string blockedResource;
  std::exception_ptr callbackException;
};

void reportError(void *data, const xmlError *error) noexcept {
  auto &state = *static_cast<ParseState *>(data);
  try {
    if (!error || error->level < XML_ERR_ERROR || !state.message.empty())
      return;
    state.message = error->message ? error->message : "XML parsing failed";
    while (!state.message.empty() &&
           (state.message.back() == '\n' || state.message.back() == '\r'))
      state.message.pop_back();
    state.location.file = error->file ? error->file : state.source.string();
    if (error->line > 0)
      state.location.line = static_cast<std::size_t>(error->line);
  } catch (...) {
    state.callbackException = std::current_exception();
  }
}

xmlParserErrors loadResource(void *data, const char *url, const char *,
                             xmlResourceType type, xmlParserInputFlags,
                             xmlParserInput **out) noexcept {
  auto &state = *static_cast<ParseState *>(data);
  *out = nullptr;
  try {
    if (!state.dtdRoot ||
        (type != XML_RESOURCE_DTD && type != XML_RESOURCE_PARAMETER_ENTITY)) {
      state.blockedResource = url ? url : "<unknown>";
      return XML_IO_EACCES;
    }
    const auto path = std::filesystem::canonical(localPath(url));
    const auto relative = path.lexically_relative(*state.dtdRoot);
    if (relative.empty() || relative.is_absolute() ||
        *relative.begin() == ".." || !std::filesystem::is_regular_file(path)) {
      state.blockedResource = url;
      return XML_IO_EACCES;
    }
    const auto contents = readFile(path);
    const auto uri = fileUri(path);
    *out = xmlNewInputFromMemory(uri.c_str(), contents.data(), contents.size(),
                                 static_cast<xmlParserInputFlags>(0));
    return *out ? XML_ERR_OK : XML_ERR_NO_MEMORY;
  } catch (...) {
    state.callbackException = std::current_exception();
    return XML_IO_EACCES;
  }
}

std::string checkedName(std::string_view value) {
  std::string name(value);
  if (name.empty() || name.find('\0') != std::string::npos ||
      xmlValidateNCName(reinterpret_cast<const xmlChar *>(name.c_str()), 0) != 0 ||
      name == "xmlns")
    throw Error("Invalid or unsupported XML name: " + name);
  return name;
}

std::string checkedText(std::string_view value) {
  if (value.size() > INT_MAX) throw Error("XML text exceeds backend capacity");
  std::string text(value);
  const auto *bytes = reinterpret_cast<const unsigned char *>(text.c_str());
  if (text.find('\0') != std::string::npos || !xmlCheckUTF8(bytes))
    throw Error("Invalid UTF-8 XML text");
  std::size_t offset = 0;
  while (offset < text.size()) {
    int length = static_cast<int>(text.size() - offset);
    const int code = xmlGetUTF8Char(bytes + offset, &length);
    if (length <= 0 || !(code == 9 || code == 10 || code == 13 ||
        (code >= 0x20 && code <= 0xD7FF) ||
        (code >= 0xE000 && code <= 0xFFFD) ||
        (code >= 0x10000 && code <= 0x10FFFF)))
      throw Error("Invalid XML character");
    offset += static_cast<std::size_t>(length);
  }
  return text;
}

class Libxml2Storage final : public detail::WritableStorage {
public:
  Libxml2Storage(NativeDocument document, std::filesystem::path source)
      : document_(std::move(document)), source_(std::move(source)) {
    std::vector<xmlNode *> pending{xmlDocGetRootElement(document_.get())};
    while (!pending.empty()) {
      auto *node = pending.back();
      pending.pop_back();
      ids_.emplace(node, nodes_.size());
      nodes_.push_back(node);
      for (auto *child = node->last; child; child = child->prev)
        if (child->type == XML_ELEMENT_NODE)
          pending.push_back(child);
    }
  }

  detail::NodeId root() const override { return 0; }

  std::string name(detail::NodeId id) const override {
    const auto *node = nodes_.at(id);
    return qualifiedName(node->name, node->ns);
  }

  std::optional<std::string> attribute(detail::NodeId id,
                                       std::string_view name) const override {
    const auto *node = nodes_.at(id);
    for (auto *attr = node->properties; attr; attr = attr->next)
      if (qualifiedName(attr->name, attr->ns) == name) {
        XmlString value(xmlNodeListGetString(document_.get(), attr->children, 1));
        return string(value.get());
      }
    for (auto *ns = node->nsDef; ns; ns = ns->next) {
      const auto declaration = ns->prefix ? "xmlns:" + string(ns->prefix) : "xmlns";
      if (declaration == name)
        return string(ns->href);
    }
    return std::nullopt;
  }

  std::string text(detail::NodeId id) const override {
    std::string result;
    for (auto *child = nodes_.at(id)->children; child; child = child->next)
      if (child->type == XML_TEXT_NODE || child->type == XML_CDATA_SECTION_NODE)
        result += string(child->content);
    return result;
  }

  std::vector<detail::NodeId> children(detail::NodeId id) const override {
    std::vector<detail::NodeId> result;
    for (auto *child = nodes_.at(id)->children; child; child = child->next)
      if (child->type == XML_ELEMENT_NODE)
        result.push_back(ids_.at(child));
    return result;
  }

  std::optional<detail::NodeId> parent(detail::NodeId id) const override {
    const auto found = ids_.find(nodes_.at(id)->parent);
    return found == ids_.end() ? std::nullopt
                               : std::optional<detail::NodeId>(found->second);
  }

  std::optional<detail::NodeId> firstChild(detail::NodeId id) const override {
    for (auto *child = nodes_.at(id)->children; child; child = child->next)
      if (child->type == XML_ELEMENT_NODE)
        return ids_.at(child);
    return std::nullopt;
  }

  std::optional<detail::NodeId> nextSibling(detail::NodeId id) const override {
    for (auto *sibling = nodes_.at(id)->next; sibling; sibling = sibling->next)
      if (sibling->type == XML_ELEMENT_NODE)
        return ids_.at(sibling);
    return std::nullopt;
  }

  SourceLocation location(detail::NodeId id) const override {
    auto *node = nodes_.at(id);
    XmlString path(xmlGetNodePath(node));
    SourceLocation result{source_.string(), string(path.get()), {}};
    const auto line = xmlGetLineNo(node);
    if (line > 0) result.line = static_cast<std::size_t>(line);
    return result;
  }

  detail::NodeId appendChild(detail::NodeId id, std::string_view name) override {
    const auto valid = checkedName(name);
    auto *node = xmlNewDocNode(document_.get(), nullptr,
        reinterpret_cast<const xmlChar *>(valid.c_str()), nullptr);
    if (!node) throw Error("Cannot allocate XML element", location(id));
    const auto newId = nodes_.size();
    try {
      nodes_.push_back(node);
      try { ids_.emplace(node, newId); }
      catch (...) { nodes_.pop_back(); throw; }
    } catch (...) { xmlFreeNode(node); throw; }
    xmlAddChild(nodes_.at(id), node);
    return newId;
  }

  void setAttribute(detail::NodeId id, std::string_view name,
                    std::string_view value) override {
    const auto valid = checkedName(name);
    const auto text = checkedText(value);
    if (!xmlSetProp(nodes_.at(id),
                    reinterpret_cast<const xmlChar *>(valid.c_str()),
                    reinterpret_cast<const xmlChar *>(text.c_str())))
      throw Error("Cannot set XML attribute", location(id));
  }

  void appendText(detail::NodeId id, std::string_view value) override {
    const auto text = checkedText(value);
    auto *node = xmlNewDocTextLen(document_.get(),
        reinterpret_cast<const xmlChar *>(text.c_str()), static_cast<int>(text.size()));
    if (!node) throw Error("Cannot allocate XML text", location(id));
    xmlAddChild(nodes_.at(id), node);
  }

  std::string serialize(bool pretty) const override {
    xmlChar *buffer = nullptr;
    int size = 0;
    xmlDocDumpFormatMemoryEnc(document_.get(), &buffer, &size, "UTF-8", pretty ? 1 : 0);
    XmlString owned(buffer);
    if (!owned || size <= 0) throw Error("Cannot serialize XML document");
    return std::string(reinterpret_cast<const char *>(owned.get()),
                       static_cast<std::size_t>(size));
  }

private:
  NativeDocument document_;
  std::filesystem::path source_;
  std::vector<xmlNode *> nodes_;
  std::unordered_map<const xmlNode *, detail::NodeId> ids_;
};

} // namespace

WritableDocument Libxml2Backend::createDocument(std::string_view rootName) const {
  const auto name = checkedName(rootName);
  NativeDocument document(xmlNewDoc(BAD_CAST "1.0"), xmlFreeDoc);
  if (!document) throw Error("Cannot allocate XML document");
  auto *root = xmlNewDocNode(document.get(), nullptr,
      reinterpret_cast<const xmlChar *>(name.c_str()), nullptr);
  if (!root) throw Error("Cannot allocate XML root");
  xmlDocSetRootElement(document.get(), root);
  return WritableDocument(std::make_shared<Libxml2Storage>(std::move(document),
                                                           std::filesystem::path{}));
}

Document Libxml2Backend::parseFile(const std::filesystem::path &file,
                                   const ParseOptions &options) const {
  return parseString(readFile(file), file, options);
}

Document Libxml2Backend::parseString(std::string_view source,
                                     const std::filesystem::path &sourceFile,
                                     const ParseOptions &options) const {
  const auto path =
      std::filesystem::absolute(sourceFile.empty() ? "memory.xml" : sourceFile)
          .lexically_normal();
  if (source.size() > INT_MAX)
    throw Error("XML input exceeds parser capacity", {path.string(), {}, {}});
  ParseState state;
  state.source = path;
  if (options.dtd == DtdPolicy::Disabled && options.dtdDirectory)
    throw Error("dtdDirectory requires DTD loading or validation",
                {path.string(), {}, {}});
  if (options.dtd == DtdPolicy::Validate && !xmlHasFeature(XML_WITH_VALID))
    throw Error("libxml2 was built without DTD validation");
  if (options.dtd != DtdPolicy::Disabled && options.dtdDirectory) {
    std::error_code error;
    auto root = std::filesystem::canonical(*options.dtdDirectory, error);
    if (error || !std::filesystem::is_directory(root))
      throw Error("DTD directory does not exist",
                  {options.dtdDirectory->string(), {}, {}});
    state.dtdRoot = std::move(root);
  }
  NativeContext context(xmlNewParserCtxt(), xmlFreeParserCtxt);
  if (!context)
    throw Error("Cannot allocate XML parser context", {path.string(), {}, {}});
  xmlCtxtSetErrorHandler(context.get(), reportError, &state);
  xmlCtxtSetResourceLoader(context.get(), loadResource, &state);
  int flags = XML_PARSE_NONET | XML_PARSE_NOENT | XML_PARSE_BIG_LINES;
  if (options.dtd != DtdPolicy::Disabled)
    flags |= XML_PARSE_DTDLOAD | XML_PARSE_DTDATTR;
  if (options.dtd == DtdPolicy::Validate)
    flags |= XML_PARSE_DTDVALID;
  const auto uri = fileUri(path);
  NativeDocument document(xmlCtxtReadMemory(context.get(), source.data(),
                                            static_cast<int>(source.size()),
                                            uri.c_str(), nullptr, flags),
                          xmlFreeDoc);
  if (state.callbackException) {
    try {
      std::rethrow_exception(state.callbackException);
    } catch (const std::exception &error) {
      throw Error(std::string("XML resource load failed: ") + error.what(),
                  {path.string(), {}, {}});
    } catch (...) {
      throw Error("XML parser callback failed", {path.string(), {}, {}});
    }
  }
  if (!state.blockedResource.empty())
    throw Error("XML resource is not permitted: " + state.blockedResource,
                {path.string(), {}, {}});
  if (!document || !state.message.empty() ||
      xmlCtxtGetStatus(context.get()) != 0 ||
      !xmlDocGetRootElement(document.get()))
    throw Error(state.message.empty() ? "XML parsing failed" : state.message,
                state.location.file.empty()
                    ? SourceLocation{path.string(), {}, {}}
                    : state.location);
  return Document(std::make_shared<Libxml2Storage>(std::move(document), path));
}

} // namespace grains::xml
