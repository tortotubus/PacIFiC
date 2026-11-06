// Parse CCCL before the legacy common headers import names from std.
#include <cub/device/device_select.cuh>
#include <thrust/iterator/counting_iterator.h>

#include "CollisionCompactionCUDA.hh"
#include "GrainsCUDAUtils.hh"

namespace grains::cuda::detail
{
    std::size_t collisionCompactionStorageBytes(unsigned int count)
    {
        if(count == 0)
            return 0;
        std::size_t                             bytes = 0;
        thrust::counting_iterator<unsigned int> input(0);
        cudaErrCheck(cub::DeviceSelect::Flagged(nullptr,
                                                bytes,
                                                input,
                                                static_cast<const std::uint8_t*>(nullptr),
                                                static_cast<unsigned int*>(nullptr),
                                                static_cast<int*>(nullptr),
                                                static_cast<int>(count)));
        return bytes;
    }
    unsigned int compactCollisionFlags(const std::uint8_t* flags,
                                       unsigned int        count,
                                       unsigned int*       indices,
                                       int*                countDevice,
                                       const int*          countHost,
                                       void*               scratch,
                                       std::size_t         scratchBytes)
    {
        if(count == 0)
            return 0;
        // Stable selection of counting indices preserves original pair IDs.
        thrust::counting_iterator<unsigned int> input(0);
        cudaErrCheck(cub::DeviceSelect::Flagged(scratch,
                                                scratchBytes,
                                                input,
                                                flags,
                                                indices,
                                                countDevice,
                                                static_cast<int>(count)));
        cudaErrCheck(cudaDeviceSynchronize());
        return static_cast<unsigned int>(*countHost);
    }
}
