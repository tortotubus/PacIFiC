// Include CCCL before the legacy common headers, which import std names.
#include <cub/device/device_select.cuh>
#include <thrust/iterator/counting_iterator.h>

#include "ForceCompactionCUDA.hh"
#include "ContactInfo.hh"
#include "GrainsCUDAUtils.hh"

namespace grains::cuda::detail
{
    template <typename T>
    struct NegativeOverlap
    {
        const ContactInfo<T>*    contacts;
        __host__ __device__ bool operator()(unsigned int i) const
        {
            return contacts[i].getOverlapDistance() < T(0);
        }
    };

    template <typename T>
    std::size_t compactionStorageBytes(unsigned int count)
    {
        if(count == 0)
            return 0;
        std::size_t bytes = 0;
        // A counting iterator preserves original pair IDs. Thrust's iterator also
        // works with CCCL versions that no longer expose cub::CountingInputIterator.
        thrust::counting_iterator<unsigned int> indices(0);
        cudaErrCheck(cub::DeviceSelect::If(nullptr,
                                           bytes,
                                           indices,
                                           static_cast<unsigned int*>(nullptr),
                                           static_cast<unsigned int*>(nullptr),
                                           static_cast<int>(count),
                                           NegativeOverlap<T>{nullptr}));
        return bytes;
    }

    template <typename T>
    unsigned int compactContacts(const ContactInfo<T>* contacts,
                                 unsigned int          count,
                                 unsigned int*         indices,
                                 unsigned int*         countDevice,
                                 const unsigned int*   countHost,
                                 void*                 scratch,
                                 std::size_t           scratchBytes)
    {
        if(count == 0)
            return 0;
        thrust::counting_iterator<unsigned int> input(0);
        cudaErrCheck(cub::DeviceSelect::If(scratch,
                                           scratchBytes,
                                           input,
                                           indices,
                                           countDevice,
                                           static_cast<int>(count),
                                           NegativeOverlap<T>{contacts}));
        cudaErrCheck(cudaDeviceSynchronize());
        return *countHost;
    }

    template std::size_t  compactionStorageBytes<float>(unsigned int);
    template std::size_t  compactionStorageBytes<double>(unsigned int);
    template unsigned int compactContacts<float>(const ContactInfo<float>*,
                                                 unsigned int,
                                                 unsigned int*,
                                                 unsigned int*,
                                                 const unsigned int*,
                                                 void*,
                                                 std::size_t);
    template unsigned int compactContacts<double>(const ContactInfo<double>*,
                                                  unsigned int,
                                                  unsigned int*,
                                                  unsigned int*,
                                                  const unsigned int*,
                                                  void*,
                                                  std::size_t);
}
