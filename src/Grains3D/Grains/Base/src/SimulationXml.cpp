#include "SimulationXml.hh"
#include "Libxml2Backend.hh"
#include "Data.hh"

namespace grains::xml {
Document readSimulationFile(const std::filesystem::path& file) {
    ParseOptions options;
    options.dtd = DtdPolicy::LoadDefaults;
    options.dtdDirectory = get_dtd_dir();
    return Libxml2Backend{}.parseFile(file, options);
}
}
