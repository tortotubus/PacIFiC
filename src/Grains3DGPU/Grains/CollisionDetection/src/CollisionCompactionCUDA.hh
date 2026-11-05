#ifndef GRAINS_COLLISION_COMPACTION_CUDA_HH
#define GRAINS_COLLISION_COMPACTION_CUDA_HH

#include <cstddef>
#include <cstdint>
namespace grains::cuda::detail
{
    std::size_t  collisionCompactionStorageBytes(unsigned int count);
    unsigned int compactCollisionFlags(const std::uint8_t* flags,
                                       unsigned int        count,
                                       unsigned int*       indices,
                                       int*                countDevice,
                                       const int*          countHost,
                                       void*               scratch,
                                       std::size_t         scratchBytes);
}
#endif
