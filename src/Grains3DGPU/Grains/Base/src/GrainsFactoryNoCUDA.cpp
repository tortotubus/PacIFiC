#include "GrainsFactoryBackend.hh"

#include <stdexcept>

namespace grains::detail
{
template <typename T>
Grains<T>* createAcceleratedGrains(const std::string&)
{
    throw std::runtime_error("This Grains3DGPU build has no CUDA backend; use Type=Standard "
                             "or rebuild with PACIFIC_USE_CUDA=ON.");
}

template Grains<float>* createAcceleratedGrains(const std::string&);
template Grains<double>* createAcceleratedGrains(const std::string&);
}  // namespace grains::detail
