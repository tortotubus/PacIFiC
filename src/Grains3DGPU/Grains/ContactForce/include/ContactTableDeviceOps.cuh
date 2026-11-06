#ifndef CONTACT_TABLE_DEVICE_OPS_CUH
#define CONTACT_TABLE_DEVICE_OPS_CUH

#include "ContactTable.hh"

namespace grains::cuda
{
struct ContactTableOps
{
    template <typename T>
    __DEVICE__ static bool find(const ContactMemoryView<T>& view, IndexPair key, uint& index)
    {
        if(view.m_capacity == 0 || view.m_table == nullptr)
            return false;

        const uint hash = primeHash(key) % view.m_capacity;
        for(uint i = 0; i < view.m_capacity; ++i)
        {
            const uint          slot  = (hash + i) % view.m_capacity;
            const ContactEntry& entry = view.m_table[slot];
            if(entry.m_valid == 0)
                return false;
            if((entry.m_valid == 1 || entry.m_valid == 2) && entry.m_key.x == key.x
               && entry.m_key.y == key.y)
            {
                __threadfence();
                index = entry.m_index;
                return true;
            }
        }
        return false;
    }

    template <typename T>
    __DEVICE__ static bool findOrInsert(ContactMemoryView<T> view, IndexPair key, uint& index)
    {
        if(view.m_capacity == 0 || view.m_table == nullptr || view.m_nextIndex == nullptr)
            return false;

        const uint hash = primeHash(key) % view.m_capacity;
        for(uint i = 0; i < view.m_capacity; ++i)
        {
            const uint    slot  = (hash + i) % view.m_capacity;
            ContactEntry& entry = view.m_table[slot];
            if(entry.m_valid == 0)
            {
                if(atomicCAS(&entry.m_valid, 0u, 2u) == 0u)
                {
                    entry.m_key   = key;
                    entry.m_index = atomicAdd(view.m_nextIndex, 1u);
                    __threadfence();
                    index = entry.m_index;
                    return true;
                }
            }
            else if(entry.m_valid == 1 && entry.m_key.x == key.x && entry.m_key.y == key.y)
            {
                atomicExch(&entry.m_valid, 2u);
                __threadfence();
                index = entry.m_index;
                return true;
            }
            else if(entry.m_valid == 2 && entry.m_key.x == key.x && entry.m_key.y == key.y)
            {
                __threadfence();
                index = entry.m_index;
                return true;
            }
        }
        return false;
    }

    template <typename T>
    __DEVICE__ static bool remove(ContactMemoryView<T> view, IndexPair key)
    {
        if(view.m_capacity == 0 || view.m_table == nullptr)
            return false;

        const uint hash = primeHash(key) % view.m_capacity;
        for(uint i = 0; i < view.m_capacity; ++i)
        {
            const uint    slot  = (hash + i) % view.m_capacity;
            ContactEntry& entry = view.m_table[slot];
            if(entry.m_valid == 0)
                return false;
            if(entry.m_key.x == key.x && entry.m_key.y == key.y)
            {
                atomicExch(&entry.m_valid, 0u);
                return true;
            }
        }
        return false;
    }
};
}  // namespace grains::cuda

#endif
