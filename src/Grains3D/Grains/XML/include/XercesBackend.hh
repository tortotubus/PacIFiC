#ifndef GRAINS_XERCES_BACKEND_HH
#define GRAINS_XERCES_BACKEND_HH

#include "Xml.hh"

namespace grains::xml {

// Xerces-C++ 3.x backend. Native types and process-runtime management stay
// private. Independent documents and parser instances may coexist.
class XercesBackend final : public Backend {
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
