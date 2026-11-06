#ifndef GRAINS_LIBXML2_BACKEND_HH
#define GRAINS_LIBXML2_BACKEND_HH

#include "Xml.hh"

namespace grains::xml {

// Stateless and reentrant. Native libxml2 types remain private to the source
// file, and each parsed document owns its context-independent tree storage.
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
