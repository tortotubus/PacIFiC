#include "ContactTable.hh"
#include "ContactTableCUDA.hh"
#include "GrainsCUDAError.hh"
#include "GrainsParameters.hh"

namespace grains::cuda
{
    template <typename T>
    __GLOBAL__ void markAndSweepKernel(ContactEntry*      table,
                                       ContactHistory<T>* history,
                                       uint               capacity,
                                       uint               maxContacts)
    {
        const uint idx = blockIdx.x * blockDim.x + threadIdx.x;
        if(idx >= capacity)
            return;

        ContactEntry& entry = table[idx];

        if(entry.m_valid == 1)
        {
            if(entry.m_index < maxContacts)
                history[entry.m_index] = ContactHistory<T>();
            entry.m_valid = 0;
        }
        else if(entry.m_valid == 2)
        {
            entry.m_valid = 1;
        }
    }

    unsigned int readContactCounter(const unsigned int* counter)
    {
        unsigned int value;
        cudaErrCheck(cudaMemcpy(&value, counter, sizeof(value), cudaMemcpyDeviceToHost));
        return value;
    }

    void zeroContactMemory(void* data, std::size_t bytes)
    {
        cudaErrCheck(cudaMemset(data, 0, bytes));
    }

    unsigned int* allocateContactCounter()
    {
        unsigned int* counter = nullptr;
        cudaErrCheck(cudaMalloc(&counter, sizeof(*counter)));
        resetContactCounter(counter);
        return counter;
    }

    void freeContactCounter(unsigned int* counter)
    {
        cudaFree(counter);
    }

    void resetContactCounter(unsigned int* counter)
    {
        zeroContactMemory(counter, sizeof(*counter));
    }

    template <typename T>
    void markAndSweepContacts(ContactEntry*      table,
                              ContactHistory<T>* history,
                              unsigned int       capacity,
                              unsigned int       maxContacts)
    {
        uint numBlocks, numThreads;
        computeOptimalThreadsAndBlocks(capacity, GrainsParameters<T>::m_GPU, numBlocks, numThreads);
        markAndSweepKernel<T><<<numBlocks, numThreads>>>(table, history, capacity, maxContacts);
        cudaErrCheck(cudaGetLastError());
    }
}  // namespace grains::cuda

#include "ContactTableImpl.hh"

template class ContactHashTable<float, MemType::DEVICE>;
template class ContactHashTable<double, MemType::DEVICE>;
