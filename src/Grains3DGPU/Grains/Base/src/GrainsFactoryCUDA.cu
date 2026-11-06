#include "GrainsFactoryBackend.hh"

#include "GrainsGPU.hh"
#include "GrainsTestDev.hh"

namespace grains::detail
{
template <typename T>
Grains<T>* createAcceleratedGrains(const std::string& option)
{
    if(option == "GPU")
        return new GrainsGPU<T>();
    return new GrainsTestDev<T>();
}

template Grains<float>* createAcceleratedGrains(const std::string&);
template Grains<double>* createAcceleratedGrains(const std::string&);
}  // namespace grains::detail
