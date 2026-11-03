#ifndef GRAINS_SIMULATION_XML_HH
#define GRAINS_SIMULATION_XML_HH

#include "Xml.hh"

namespace grains::xml {
Document readSimulationFile(const std::filesystem::path& file);
}

#endif
