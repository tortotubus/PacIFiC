// Keep CUB independent of Basic.hh and its legacy namespace imports.
#include "GrainsCUDAError.hh"
#include "ParticleSortCUDA.hh"
#include <cub/device/device_radix_sort.cuh>

namespace grains::cuda::detail
{
    std::size_t particleSortStorageBytes(unsigned int count)
    {
        if(count == 0)
            return 0;
        std::size_t bytes = 0;
        cudaErrCheck(cub::DeviceRadixSort::SortPairs(nullptr,
                                                     bytes,
                                                     static_cast<const std::uint64_t*>(nullptr),
                                                     static_cast<std::uint64_t*>(nullptr),
                                                     static_cast<const unsigned int*>(nullptr),
                                                     static_cast<unsigned int*>(nullptr),
                                                     count));
        return bytes;
    }
    void sortParticleIndices(const std::uint64_t* codes,
                             std::uint64_t*       sortedCodes,
                             const unsigned int*  indices,
                             unsigned int*        sortedIndices,
                             unsigned int         count,
                             void*                scratch,
                             std::size_t          scratchBytes)
    {
        if(count == 0)
            return;
        cudaErrCheck(cub::DeviceRadixSort::SortPairs(scratch,
                                                     scratchBytes,
                                                     codes,
                                                     sortedCodes,
                                                     indices,
                                                     sortedIndices,
                                                     count));
        cudaErrCheck(cudaDeviceSynchronize());
    }
}
