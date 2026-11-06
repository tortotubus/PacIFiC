#ifndef CONTACT_TABLE_HOST_OPS_HH
#define CONTACT_TABLE_HOST_OPS_HH

#include "ContactTable.hh"

namespace grains::cpu
{
struct ContactTableOps
{
    template <typename T>
    static bool find(const ContactMemoryView<T>& view, IndexPair key, uint& index)
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
                index = entry.m_index;
                return true;
            }
        }
        return false;
    }

    template <typename T>
    static bool findOrInsert(ContactMemoryView<T> view, IndexPair key, uint& index)
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
                entry.m_valid = 2;
                entry.m_key   = key;
                entry.m_index = (*view.m_nextIndex)++;
                index         = entry.m_index;
                return true;
            }
            if(entry.m_valid == 1 && entry.m_key.x == key.x && entry.m_key.y == key.y)
            {
                entry.m_valid = 2;
                index         = entry.m_index;
                return true;
            }
            if(entry.m_valid == 2 && entry.m_key.x == key.x && entry.m_key.y == key.y)
            {
                index = entry.m_index;
                return true;
            }
        }
        return false;
    }

    template <typename T>
    static bool remove(ContactMemoryView<T> view, IndexPair key)
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
                entry.m_valid = 0;
                return true;
            }
        }
        return false;
    }
};
}  // namespace grains::cpu

#endif
