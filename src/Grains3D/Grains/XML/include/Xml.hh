#ifndef GRAINS_XML_HH
#define GRAINS_XML_HH

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <locale>
#include <iomanip>
#include <limits>
#include <cmath>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace grains::xml {

struct SourceLocation {
    std::string file;
    std::string path;
    std::optional<std::size_t> line;
};

class Error : public std::runtime_error {
public:
    Error(std::string message, SourceLocation location = {})
        : std::runtime_error(describe(message, location)), location_(std::move(location)) {}

    const SourceLocation& location() const noexcept { return location_; }

private:
    static std::string describe(const std::string& message, const SourceLocation& where) {
        std::string prefix = where.file;
        if (where.line) prefix += ":" + std::to_string(*where.line);
        if (!where.path.empty()) prefix += " " + where.path;
        return prefix.empty() ? message : prefix + ": " + message;
    }
    SourceLocation location_;
};

namespace detail {

using NodeId = std::size_t;

// Backend extension point. IDs are document-local and stable for the lifetime
// of this tree; writable backends must retain existing element IDs. No particular ID is reserved as a null sentinel.
// Strings are UTF-8, names are qualified names, children are elements only in
// document order, and text concatenates direct text/CDATA children (untrimmed).
// Implementations own all native resources and must not expose native types.
class DocumentStorage {
public:
    virtual ~DocumentStorage() = default;
    virtual NodeId root() const = 0;
    virtual std::string name(NodeId node) const = 0;
    virtual std::optional<std::string> attribute(NodeId node, std::string_view name) const = 0;
    virtual std::string text(NodeId node) const = 0;
    virtual std::vector<NodeId> children(NodeId node) const = 0;
    virtual SourceLocation location(NodeId node) const = 0;
    // Default keeps existing storage implementations source compatible.
    // Native tree backends should override this with constant-time lookup.
    virtual std::optional<NodeId> parent(NodeId node) const {
        if (node == root()) return std::nullopt;
        std::vector<NodeId> pending{root()};
        while (!pending.empty()) {
            const auto current = pending.back();
            pending.pop_back();
            for (auto child : children(current)) {
                if (child == node) return current;
                pending.push_back(child);
            }
        }
        return std::nullopt;
    }
    virtual std::optional<NodeId> firstChild(NodeId node) const {
        const auto elements = children(node);
        if (elements.empty()) return std::nullopt;
        return elements.front();
    }
    virtual std::optional<NodeId> nextSibling(NodeId node) const {
        const auto owner = parent(node);
        if (!owner) return std::nullopt;
        const auto elements = children(*owner);
        for (std::size_t i = 0; i + 1 < elements.size(); ++i)
            if (elements[i] == node) return elements[i + 1];
        return std::nullopt;
    }
};

// Separate capability: parsed documents remain read-only. Writable handles
// support append-only trees so existing element IDs never become invalid.
class WritableStorage : public DocumentStorage {
public:
    virtual NodeId appendChild(NodeId parent, std::string_view name) = 0;
    virtual void setAttribute(NodeId node, std::string_view name, std::string_view value) = 0;
    virtual void appendText(NodeId node, std::string_view value) = 0;
    virtual std::string serialize(bool pretty) const = 0;
};

template<class T>
std::string format(T value) {
    static_assert(std::is_arithmetic_v<T>, "XML formatting supports arithmetic types");
    static_assert(!std::is_same_v<T, char> && !std::is_same_v<T, signed char> &&
                  !std::is_same_v<T, unsigned char>, "Use an integer type");
    if constexpr (std::is_floating_point_v<T>)
        if (!std::isfinite(value)) throw Error("Non-finite XML numeric value");
    std::ostringstream output;
    output.imbue(std::locale::classic());
    if constexpr (std::is_floating_point_v<T>)
        output << std::setprecision(std::numeric_limits<T>::max_digits10);
    if constexpr (std::is_same_v<T, bool>) output << std::boolalpha;
    output << value;
    return output.str();
}

template<class T>
T convert(const std::string& value, const SourceLocation& where) {
    static_assert(std::is_same_v<T, std::string> || std::is_arithmetic_v<T>,
                  "XML conversion supports strings and arithmetic types");
    if constexpr (std::is_same_v<T, std::string>) {
        return value;
    } else {
        std::istringstream input(value);
        input.imbue(std::locale::classic());
        input >> std::ws;
        if constexpr (std::is_same_v<T, bool>) {
            std::string token;
            input >> token;
            input >> std::ws;
            if (input.eof()) {
                if (token == "true" || token == "True" || token == "1") return true;
                if (token == "false" || token == "False" || token == "0") return false;
            }
        } else {
            // Stream extraction into unsigned integers otherwise accepts -1.
            if constexpr (std::is_unsigned_v<T>) {
                if (input.peek() == '-') throw Error("Negative unsigned value '" + value + "'", where);
            }
            T result{};
            // Character types must be parsed as numbers, not single letters.
            if constexpr (std::is_same_v<T, char> || std::is_same_v<T, signed char> ||
                          std::is_same_v<T, unsigned char>) {
                static_assert(!std::is_same_v<T, T>, "Use an integer type rather than a character type");
            }
            if (input >> result) {
                input >> std::ws;
                if (input.eof()) return result;
            }
        }
        throw Error("Invalid value '" + value + "'", where);
    }
}

} // namespace detail

class Document;

// A non-null, read-only handle retaining ownership of its document storage.
// A view of writable storage observes subsequent writes.
class Element {
public:
    std::string name() const { return storage_->name(id_); }
    SourceLocation location() const { return storage_->location(id_); }
    std::optional<std::string> attribute(std::string_view name) const {
        return storage_->attribute(id_, name);
    }
    bool hasAttribute(std::string_view name) const { return attribute(name).has_value(); }
    std::string text() const { return storage_->text(id_); }

    std::vector<Element> children() const {
        std::vector<Element> result;
        for (auto id : storage_->children(id_)) result.push_back(Element(storage_, id));
        return result;
    }
    std::vector<Element> children(std::string_view name) const {
        std::vector<Element> result;
        for (const auto& element : children())
            if (element.name() == name) result.push_back(element);
        return result;
    }
    // First matching direct child. Use requireChild for singleton sections.
    std::optional<Element> child(std::string_view name) const {
        for (const auto& element : children())
            if (element.name() == name) return element;
        return std::nullopt;
    }
    Element requireChild(std::string_view name) const {
        auto matches = children(name);
        if (matches.size() != 1)
            throw Error("Expected exactly one child '" + std::string(name) + "'", location());
        return matches.front();
    }
    std::optional<Element> firstChild() const {
        auto id = storage_->firstChild(id_);
        if (!id) return std::nullopt;
        return Element(storage_, *id);
    }
    std::optional<Element> parent() const {
        auto id = storage_->parent(id_);
        if (!id) return std::nullopt;
        return Element(storage_, *id);
    }
    std::optional<Element> nextSibling() const {
        auto id = storage_->nextSibling(id_);
        if (!id) return std::nullopt;
        return Element(storage_, *id);
    }
    // Descendants exclude this element and retain document order. This is
    // deliberately separate from the direct-child lookup above.
    std::vector<Element> descendants() const {
        std::vector<Element> result;
        auto pending = children();
        std::reverse(pending.begin(), pending.end());
        while (!pending.empty()) {
            auto element = pending.back();
            pending.pop_back();
            result.push_back(element);
            auto nested = element.children();
            pending.insert(pending.end(), nested.rbegin(), nested.rend());
        }
        return result;
    }
    std::vector<Element> descendants(std::string_view name) const {
        std::vector<Element> result;
        for (const auto& element : descendants())
            if (element.name() == name) result.push_back(element);
        return result;
    }
    std::optional<Element> descendant(std::string_view name) const {
        for (const auto& element : descendants())
            if (element.name() == name) return element;
        return std::nullopt;
    }
    // Explicit search-scope names used by migrated simulation readers.
    std::optional<Element> getChild(std::string_view name) const { return child(name); }
    std::optional<Element> getDescendant(std::string_view name) const { return descendant(name); }
    std::vector<Element> getChildren() const { return children(); }
    std::vector<Element> getChildren(std::string_view name) const { return children(name); }
    std::vector<Element> getDescendants() const { return descendants(); }
    std::vector<Element> getDescendants(std::string_view name) const { return descendants(name); }
    template<class T>
    std::optional<T> attributeAs(std::string_view name) const {
        auto value = attribute(name);
        if (!value) return std::nullopt;
        return convertAttribute<T>(name, *value);
    }
    template<class T>
    T requireAttribute(std::string_view name) const {
        auto value = attribute(name);
        if (!value) throw Error("Missing attribute '" + std::string(name) + "'", location());
        return convertAttribute<T>(name, *value);
    }
    template<class T>
    T attributeOr(std::string_view name, T fallback) const {
        auto value = attribute(name);
        return value ? convertAttribute<T>(name, *value) : std::move(fallback);
    }
    template<class T>
    T value() const { return detail::convert<T>(text(), location()); }

private:
    friend class Document;
    friend class WritableElement;
    Element(std::shared_ptr<const detail::DocumentStorage> storage, detail::NodeId id)
        : storage_(std::move(storage)), id_(id) {}
    template<class T>
    T convertAttribute(std::string_view name, const std::string& value) const {
        auto where = location();
        where.path += "/@" + std::string(name);
        return detail::convert<T>(value, where);
    }
    std::shared_ptr<const detail::DocumentStorage> storage_;
    detail::NodeId id_;
};

// Nullable owning handles replace the nullable DOMNode pointers used by legacy
// simulation readers. Element itself remains non-null. Required access throws
// rather than dereferencing an empty optional; presence tests stay explicit.
using OptionalElement = std::optional<Element>;

inline const Element& requireElement(const Element& element) noexcept { return element; }
inline const Element& requireElement(const OptionalElement& element) {
    if (!element) throw Error("Missing XML element");
    return *element;
}
inline bool hasAttribute(const Element& element, std::string_view name) {
    return element.hasAttribute(name);
}
inline bool hasAttribute(const OptionalElement& element, std::string_view name) {
    return element && element->hasAttribute(name);
}

class Document {
public:
    explicit Document(std::shared_ptr<const detail::DocumentStorage> storage)
        : storage_(std::move(storage)) {
        if (!storage_) throw std::invalid_argument("XML document storage is null");
    }
    Element root() const { return Element(storage_, storage_->root()); }
private:
    std::shared_ptr<const detail::DocumentStorage> storage_;
};

// Copyable owning handles; writes are visible through every handle to the
// same document. Concurrent access requires external synchronization.
class WritableElement {
public:
    Element read() const { return Element(storage_, id_); }
    WritableElement appendChild(std::string_view name) const {
        return WritableElement(storage_, storage_->appendChild(id_, name));
    }
    void setAttribute(std::string_view name, std::string_view value) const {
        storage_->setAttribute(id_, name, value);
    }
    template<class T, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
    void setAttribute(std::string_view name, T value) const {
        setAttribute(name, detail::format(value));
    }
    // Append literal character data, escaping XML metacharacters. Existing
    // children and text are retained, including mixed-content ordering.
    void appendText(std::string_view value) const { storage_->appendText(id_, value); }
    template<class T, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
    void appendText(T value) const { appendText(detail::format(value)); }
private:
    friend class WritableDocument;
    WritableElement(std::shared_ptr<detail::WritableStorage> storage, detail::NodeId id)
        : storage_(std::move(storage)), id_(id) {}
    std::shared_ptr<detail::WritableStorage> storage_;
    detail::NodeId id_;
};

class WritableDocument {
public:
    explicit WritableDocument(std::shared_ptr<detail::WritableStorage> storage)
        : storage_(std::move(storage)) {
        if (!storage_) throw std::invalid_argument("XML document storage is null");
    }
    WritableElement root() const { return WritableElement(storage_, storage_->root()); }
    Document read() const { return Document(storage_); }
    // UTF-8 with XML declaration. Formatting is backend-defined; mixed text
    // content must be preserved. writeFile truncates an existing file.
    std::string toString(bool pretty = true) const { return storage_->serialize(pretty); }
    void writeFile(const std::filesystem::path& file, bool pretty = true) const;
private:
    std::shared_ptr<detail::WritableStorage> storage_;
};

inline void WritableDocument::writeFile(const std::filesystem::path& file, bool pretty) const {
    const auto data = toString(pretty);
    if (data.size() > static_cast<std::size_t>(std::numeric_limits<std::streamsize>::max()))
        throw Error("XML output exceeds stream capacity", {file.string(), {}, {}});
    std::ofstream output(file, std::ios::binary | std::ios::trunc);
    if (!output) throw Error("Cannot open XML output", {file.string(), {}, {}});
    output.write(data.data(), static_cast<std::streamsize>(data.size()));
    output.close();
    if (!output) throw Error("Cannot write XML output", {file.string(), {}, {}});
}

// Disabled skips external DTD loading and validation. Internal declarations
// may still be parsed as part of well-formed XML (e.g. internal text entities).
enum class DtdPolicy { Disabled, Validate, LoadDefaults };

struct ParseOptions {
    DtdPolicy dtd = DtdPolicy::Disabled;
    // When loading/validating, DTD resources may only resolve inside this local
    // directory. With no directory, only an internal DTD may be used.
    // Network access is forbidden. Unsupported policies must throw.
    std::optional<std::filesystem::path> dtdDirectory;
};

class Backend {
public:
    virtual ~Backend() = default;
    virtual WritableDocument createDocument(std::string_view rootName) const {
        (void)rootName;
        throw Error("XML writing is unsupported by this backend");
    }
    // Throw Error on parse failure or unsupported options. Parsing must not
    // return partially recovered documents or depend on a global parser.
    virtual Document parseFile(const std::filesystem::path& file,
                               const ParseOptions& options = {}) const = 0;
    // sourceFile supplies diagnostics and a base path for relative DTD paths.
    virtual Document parseString(std::string_view source,
                                 const std::filesystem::path& sourceFile = "memory.xml",
                                 const ParseOptions& options = {}) const {
        (void)source; (void)options;
        throw Error("In-memory parsing is unsupported by this backend", {sourceFile.string(), {}, {}});
    }
};

} // namespace grains::xml

#endif
