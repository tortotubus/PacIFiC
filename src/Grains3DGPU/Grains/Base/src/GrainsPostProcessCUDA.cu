#include "GrainsPostProcessImpl.hh"

template void Grains<float>::postProcess<MemType::DEVICE>(
    const std::unique_ptr<ComponentManager<float, MemType::DEVICE>>&);
template void Grains<double>::postProcess<MemType::DEVICE>(
    const std::unique_ptr<ComponentManager<double, MemType::DEVICE>>&);
