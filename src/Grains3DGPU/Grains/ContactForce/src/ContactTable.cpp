#include "ContactTable.hh"

// -------------------------------------------------------------------------------------------------
// CUDA kernel for mark-and-sweep cleanup
template <typename T>
__GLOBAL__ void markAndSweepKernel(ContactEntry* table, uint capacity)
{
    uint idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx >= capacity)
        return;

    if(table[idx].m_valid == ContactEntry::ACTIVE)
        table[idx].m_valid = ContactEntry::TOMBSTONE;
}

// -------------------------------------------------------------------------------------------------
// Default constructor
template <typename T, MemType M>
ContactHashTable<T, M>::ContactHashTable()
    : m_historyData()
    , m_table()
    , m_capacity(0)
    , m_maxContacts(0)
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
// Gets a complete view for passing to kernels
template <typename T, MemType M>
ContactMemoryView<T> ContactHashTable<T, M>::getView()
{
    return ContactMemoryView<T>(m_historyData.getData(), m_table.getData(), m_capacity);
}

// -------------------------------------------------------------------------------------------------
// Counts the number of ACTIVE entries (host-side, O(capacity))
template <typename T, MemType M>
uint ContactHashTable<T, M>::countActive() const
{
    const ContactEntry* table = m_table.getData();
    uint                count = 0;
    for(uint i = 0; i < m_capacity; ++i)
        if(table[i].m_valid == ContactEntry::ACTIVE)
            ++count;
    return count;
}

// -------------------------------------------------------------------------------------------------
// Allocates memory for both hash table and history data
template <typename T, MemType M>
void ContactHashTable<T, M>::allocate(uint hashCapacity, uint maxContacts)
{
    // Deallocate existing memory if any
    if(m_table.getSize() > 0 || m_historyData.getSize() > 0)
        deallocate();

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
        cudaErrCheck(cudaMemset(m_table.getData(), 0, hashCapacity * sizeof(ContactEntry)));
    }

    // Allocate and initialize the history data (one slot per hash-table slot)
    m_historyData.initialize(hashCapacity);
    if constexpr(M == MemType::HOST)
    {
        m_historyData.fill(ContactHistory<T>());
    }
    else if constexpr(M == MemType::DEVICE)
    {
        cudaErrCheck(
            cudaMemset(m_historyData.getData(), 0, hashCapacity * sizeof(ContactHistory<T>)));
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
            cudaErrCheck(cudaMemset(m_table.getData(), 0, m_capacity * sizeof(ContactEntry)));
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
            cudaErrCheck(
                cudaMemset(m_historyData.getData(), 0, m_capacity * sizeof(ContactHistory<T>)));
        }
    }
}

// -------------------------------------------------------------------------------------------------
// Grows the contact table if newMaxContacts exceeds the current capacity.
// Existing contact history is discarded (contacts restart as new).
template <typename T, MemType M>
void ContactHashTable<T, M>::grow(uint newMaxContacts)
{
    if(newMaxContacts <= m_maxContacts)
        return;

    uint newHashCapacity = static_cast<uint>(newMaxContacts / 0.7);
    allocate(newHashCapacity, newMaxContacts);
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
        ContactEntry* table = m_table.getData();
        for(uint i = 0; i < m_capacity; ++i)
        {
            if(table[i].m_valid == ContactEntry::ACTIVE)
                table[i].m_valid = ContactEntry::TOMBSTONE;
        }
    }
    else if constexpr(M == MemType::DEVICE)
    {
        uint numBlocks, numThreads;
        computeOptimalThreadsAndBlocks(m_capacity,
                                       GrainsParameters<T>::m_GPU,
                                       numBlocks,
                                       numThreads);
        markAndSweepKernel<T><<<numBlocks, numThreads>>>(m_table.getData(), m_capacity);
        cudaErrCheck(cudaGetLastError());
    }
}

// -------------------------------------------------------------------------------------------------
// Explicit template instantiation
template class ContactHashTable<float, MemType::HOST>;
template class ContactHashTable<float, MemType::DEVICE>;
template class ContactHashTable<double, MemType::HOST>;
template class ContactHashTable<double, MemType::DEVICE>;
