#ifndef GRAINS_FORCE_COMPACTION_CUDA_HH
#define GRAINS_FORCE_COMPACTION_CUDA_HH

#include <cstddef>
template <typename T>
class ContactInfo;

namespace grains::cuda::detail
{
    template <typename T>
    std::size_t compactionStorageBytes(unsigned int count);
    template <typename T>
    unsigned int compactContacts(const ContactInfo<T>* contacts,
                                 unsigned int          count,
                                 unsigned int*         indices,
                                 unsigned int*         countDevice,
                                 const unsigned int*   countHost,
                                 void*                 scratch,
                                 std::size_t           scratchBytes);
}
#endif
