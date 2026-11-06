#ifndef GRAINS_FACTORY_BACKEND_HH
#define GRAINS_FACTORY_BACKEND_HH

#include <string>

template <typename T>
class Grains;

namespace grains::detail
{
    /** Create a GPU application, or throw when this build has no CUDA backend. */
    template <typename T>
    Grains<T>* createAcceleratedGrains(const std::string& option);
}

#endif
