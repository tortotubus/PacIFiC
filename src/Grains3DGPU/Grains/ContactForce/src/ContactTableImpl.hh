#ifndef CONTACT_TABLE_IMPL_HH
#define CONTACT_TABLE_IMPL_HH

#include "ContactTable.hh"
#include "ContactTableCUDA.hh"

// -------------------------------------------------------------------------------------------------
// Default constructor
template <typename T, MemType M>
ContactHashTable<T, M>::ContactHashTable()
    : m_historyData()
    , m_table()
    , m_capacity(0)
    , m_maxContacts(0)
    , m_nextIndex(nullptr)
{
}

// -------------------------------------------------------------------------------------------------
// Constructor with specified capacities
template <typename T, MemType M>
ContactHashTable<T, M>::ContactHashTable(uint hashCapacity, uint maxContacts)
    : m_historyData()
    , m_table()
    , m_capacity(0)
    , m_maxContacts(0)
    , m_nextIndex(nullptr)
{
    allocate(hashCapacity, maxContacts);
}

// -------------------------------------------------------------------------------------------------
// Destructor
template <typename T, MemType M>
ContactHashTable<T, M>::~ContactHashTable()
{
    deallocate();
}

// -------------------------------------------------------------------------------------------------
// Gets the pointer to the history data
template <typename T, MemType M>
const ContactHistory<T>* ContactHashTable<T, M>::getHistoryData() const
{
    return m_historyData.getData();
}

// -------------------------------------------------------------------------------------------------
// Gets mutable pointer to the history data
template <typename T, MemType M>
ContactHistory<T>* ContactHashTable<T, M>::getHistoryData()
{
    return m_historyData.getData();
}

// -------------------------------------------------------------------------------------------------
// Gets the pointer to the table data
template <typename T, MemType M>
const ContactEntry* ContactHashTable<T, M>::getTable() const
{
    return m_table.getData();
}

// -------------------------------------------------------------------------------------------------
// Gets mutable pointer to the table data
template <typename T, MemType M>
ContactEntry* ContactHashTable<T, M>::getTable()
{
    return m_table.getData();
}

// -------------------------------------------------------------------------------------------------
// Gets the capacity of the hash table
template <typename T, MemType M>
uint ContactHashTable<T, M>::getCapacity() const
{
    return m_capacity;
}

// -------------------------------------------------------------------------------------------------
// Gets the maximum number of contacts
template <typename T, MemType M>
uint ContactHashTable<T, M>::getMaxContacts() const
{
    return m_maxContacts;
}

// -------------------------------------------------------------------------------------------------
// Gets the pointer to the next index counter
template <typename T, MemType M>
uint* ContactHashTable<T, M>::getNextIndexPointer()
{
    return m_nextIndex;
}

// -------------------------------------------------------------------------------------------------
// Gets the pointer to the next index counter (const)
template <typename T, MemType M>
const uint* ContactHashTable<T, M>::getNextIndexPointer() const
{
    return m_nextIndex;
}

// -------------------------------------------------------------------------------------------------
// Gets a complete view for passing to kernels
template <typename T, MemType M>
ContactMemoryView<T> ContactHashTable<T, M>::getView()
{
    return ContactMemoryView<T>(m_historyData.getData(),
                                m_table.getData(),
                                m_capacity,
                                m_nextIndex);
}

// -------------------------------------------------------------------------------------------------
// Gets the current next index value
template <typename T, MemType M>
uint ContactHashTable<T, M>::getNextIndex() const
{
    if(m_nextIndex == nullptr)
        return 0;

    if constexpr(M == MemType::HOST)
    {
        return *m_nextIndex;
    }
    else if constexpr(M == MemType::DEVICE)
    {
        return grains::cuda::readContactCounter(m_nextIndex);
    }

    return 0;
}

// -------------------------------------------------------------------------------------------------
// Allocates memory for both hash table and history data
template <typename T, MemType M>
void ContactHashTable<T, M>::allocate(uint hashCapacity, uint maxContacts)
{
    // Deallocate existing memory if any
    if(m_table.getSize() > 0 || m_historyData.getSize() > 0 || m_nextIndex != nullptr)
    {
        deallocate();
    }

    m_capacity    = hashCapacity;
    m_maxContacts = maxContacts;

    // Allocate and initialize the hash table
    m_table.initialize(hashCapacity);
    if constexpr(M == MemType::HOST)
    {
        m_table.fill(ContactEntry());
    }
    else if constexpr(M == MemType::DEVICE)
    {
        grains::cuda::zeroContactMemory(m_table.getData(),
                                        hashCapacity * sizeof(ContactEntry));
    }

    // Allocate and initialize the history data
    m_historyData.initialize(maxContacts);
    if constexpr(M == MemType::HOST)
    {
        m_historyData.fill(ContactHistory<T>());
    }
    else if constexpr(M == MemType::DEVICE)
    {
        grains::cuda::zeroContactMemory(m_historyData.getData(),
                                        maxContacts * sizeof(ContactHistory<T>));
    }

    // Allocate the index counter
    if constexpr(M == MemType::HOST)
    {
        m_nextIndex  = new uint;
        *m_nextIndex = 0;
    }
    else if constexpr(M == MemType::DEVICE)
    {
        m_nextIndex = grains::cuda::allocateContactCounter();
    }
}

// -------------------------------------------------------------------------------------------------
// Frees memory used by hash table and history data
template <typename T, MemType M>
void ContactHashTable<T, M>::deallocate()
{
    m_table.free();
    m_historyData.free();
    m_capacity    = 0;
    m_maxContacts = 0;

    if(m_nextIndex != nullptr)
    {
        if constexpr(M == MemType::HOST)
        {
            delete m_nextIndex;
        }
        else if constexpr(M == MemType::DEVICE)
        {
            grains::cuda::freeContactCounter(m_nextIndex);
        }

        m_nextIndex = nullptr;
    }
}

// -------------------------------------------------------------------------------------------------
// Clears all entries in the hash table and resets history data
template <typename T, MemType M>
void ContactHashTable<T, M>::clear()
{
    if(m_table.getSize() > 0)
    {
        if constexpr(M == MemType::HOST)
        {
            m_table.fill(ContactEntry());
        }
        else if constexpr(M == MemType::DEVICE)
        {
            grains::cuda::zeroContactMemory(m_table.getData(),
                                            m_capacity * sizeof(ContactEntry));
        }
    }

    if(m_historyData.getSize() > 0)
    {
        if constexpr(M == MemType::HOST)
        {
            m_historyData.fill(ContactHistory<T>());
        }
        else if constexpr(M == MemType::DEVICE)
        {
            grains::cuda::zeroContactMemory(m_historyData.getData(),
                                            m_maxContacts * sizeof(ContactHistory<T>));
        }
    }

    if(m_nextIndex != nullptr)
    {
        if constexpr(M == MemType::HOST)
        {
            *m_nextIndex = 0;
        }
        else if constexpr(M == MemType::DEVICE)
        {
            grains::cuda::resetContactCounter(m_nextIndex);
        }
    }
}

// -------------------------------------------------------------------------------------------------
// Resets the index counter
template <typename T, MemType M>
void ContactHashTable<T, M>::resetIndexCounter()
{
    if(m_nextIndex != nullptr)
    {
        if constexpr(M == MemType::HOST)
        {
            *m_nextIndex = 0;
        }
        else if constexpr(M == MemType::DEVICE)
        {
            grains::cuda::resetContactCounter(m_nextIndex);
        }
    }
}

// -------------------------------------------------------------------------------------------------
// Performs mark-and-sweep cleanup
template <typename T, MemType M>
void ContactHashTable<T, M>::markAndSweep()
{
    if(m_table.getSize() == 0)
        return;

    if constexpr(M == MemType::HOST)
    {
        // CPU version
        ContactEntry*      table   = m_table.getData();
        ContactHistory<T>* history = m_historyData.getData();

        for(uint i = 0; i < m_capacity; ++i)
        {
            ContactEntry& entry = table[i];

            if(entry.m_valid == 1)
            {
                // Stale entry - remove it and reset history
                if(entry.m_index < m_maxContacts)
                {
                    history[entry.m_index] = ContactHistory<T>();
                }
                entry.m_valid = 0;
            }
            else if(entry.m_valid == 2)
            {
                // Active entry - demote to stale
                entry.m_valid = 1;
            }
        }
    }
    else if constexpr(M == MemType::DEVICE)
    {
        grains::cuda::markAndSweepContacts(m_table.getData(),
                                           m_historyData.getData(),
                                           m_capacity,
                                           m_maxContacts);
    }
}

#endif
