#ifndef GRAINS_LIBXML2_BACKEND_HH
#define GRAINS_LIBXML2_BACKEND_HH

#include "Xml.hh"

namespace grains::xml {

// Stateless reader/writer. Each parse owns its context, diagnostics, resource policy,
// and document. Parsed documents expose only read handles. Public headers never expose libxml2 types.
// Requires libxml2 >= 2.14 for per-context resource loading.
class Libxml2Backend final : public Backend {
public:
    WritableDocument createDocument(std::string_view rootName) const override;
    Document parseFile(const std::filesystem::path& file,
                       const ParseOptions& options = {}) const override;
    Document parseString(std::string_view source,
                         const std::filesystem::path& sourceFile = "memory.xml",
                         const ParseOptions& options = {}) const override;
};

} // namespace grains::xml
#endif
