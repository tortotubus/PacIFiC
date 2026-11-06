#ifndef COLLISION_DETECTION_MEMORY_BACKEND_HH
#define COLLISION_DETECTION_MEMORY_BACKEND_HH

#include <cstddef>

namespace grains::cpu
{
    std::size_t getAvailableCollisionMemory();
}

namespace grains::cuda
{
    std::size_t getAvailableCollisionMemory();
}

#endif
