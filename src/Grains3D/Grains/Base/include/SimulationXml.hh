#ifndef GRAINS_SIMULATION_XML_HH
#define GRAINS_SIMULATION_XML_HH

#include "Xml.hh"

namespace grains::xml {
// Input bootstrap uses the installed/in-tree Grains DTD directory. Load
// defaults without enabling strict validation of legacy simulation inputs.
Document readSimulationFile(const std::filesystem::path& file);
}

#endif
