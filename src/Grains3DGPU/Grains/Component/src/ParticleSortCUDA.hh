#ifndef GRAINS_PARTICLE_SORT_CUDA_HH
#define GRAINS_PARTICLE_SORT_CUDA_HH

#include <cstddef>
#include <cstdint>

namespace grains::cuda::detail
{
    std::size_t particleSortStorageBytes(unsigned int count);
    void        sortParticleIndices(const std::uint64_t* codes,
                                    std::uint64_t*       sortedCodes,
                                    const unsigned int*  indices,
                                    unsigned int*        sortedIndices,
                                    unsigned int         count,
                                    void*                scratch,
                                    std::size_t          scratchBytes);
}
#endif
