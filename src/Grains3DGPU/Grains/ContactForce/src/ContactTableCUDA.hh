#ifndef CONTACT_TABLE_CUDA_HH
#define CONTACT_TABLE_CUDA_HH

#include <cstddef>

struct ContactEntry;
template <typename T>
struct ContactHistory;

namespace grains::cuda
{
unsigned int readContactCounter(const unsigned int* counter);
void         zeroContactMemory(void* data, std::size_t bytes);
unsigned int* allocateContactCounter();
void          freeContactCounter(unsigned int* counter);
void          resetContactCounter(unsigned int* counter);

template <typename T>
void markAndSweepContacts(ContactEntry*      table,
                          ContactHistory<T>* history,
                          unsigned int       capacity,
                          unsigned int       maxContacts);
}  // namespace grains::cuda

#endif
