#ifndef _CONTACTTABLE_HH_
#define _CONTACTTABLE_HH_

#include "Basic.hh"
#include "GrainsMemBuffer.hh"
#include "GrainsUtils.hh"
#include "Vector3.hh"

// =================================================================================================
/** @brief Hash table entry for contact tracking.

    Each entry stores a key (pair of component IDs) and a state flag. The slot position in the
    table doubles as the index into the parallel ContactHistory array, so no separate index field
    is needed.

    @author A.Yazdani - 2026 - Construction */
// =================================================================================================
struct ContactEntry
{
    enum : uint
    {
        EMPTY     = 0u,
        ACTIVE    = 1u,
        TOMBSTONE = 2u
    };

    /** @name Parameters */
    //@{
    /** \brief Contact pair key (i,j) where i < j */
    uint2 m_key;
    /** \brief Entry state: 0=empty, 1=active, 2=tombstone */
    uint m_valid;
    //@}

    /** @name Constructors */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Default constructor */
    __HOSTDEVICE__
    ContactEntry()
        : m_key(make_uint2(0, 0))
        , m_valid(EMPTY)
    {
    }
    //@}
};

// =================================================================================================
/** @brief Contact history data for memory-enabled force models.

    Stores the cumulative tangential and rolling displacements for a contact pair.
    This data is accessed via the index returned by the hash table.

    @author A.Yazdani - 2026 - Construction */
// =================================================================================================
template <typename T>
struct ContactHistory
{
    /** @name Parameters */
    //@{
    /** \brief Cumulative tangential displacement (kt * delta) */
    Vector3<T> m_tangentialDisplacement;
    /** \brief Cumulative rolling friction spring-torque */
    Vector3<T> m_rollingDisplacement;
    /** \brief Previous contact normal (for plane rotation between timesteps) */
    Vector3<T> m_previousNormal;
    //@}

    /** @name Constructors */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Default constructor */
    __HOSTDEVICE__
    ContactHistory()
        : m_tangentialDisplacement(Vector3<T>(T(0), T(0), T(0)))
        , m_rollingDisplacement(Vector3<T>(T(0), T(0), T(0)))
        , m_previousNormal(Vector3<T>(T(0), T(0), T(0)))
    {
    }
    //@}
};

// =================================================================================================
/** @brief Complete contact memory view for passing to force models.

    Contains the hash table and the parallel history array. The slot position i in m_table is
    also the index into m_historyData[i], so no separate index field is needed. Contains only
    raw pointers and primitive types, making it safe to pass to kernels by value.

    @author A.Yazdani - 2026 - Construction */
// =================================================================================================
template <typename T>
struct ContactMemoryView
{
    /** @name Parameters */
    //@{
    /** \brief Pointer to flat array of contact history data (one entry per hash-table slot) */
    ContactHistory<T>* m_historyData;
    /** \brief Pointer to hash table entries */
    ContactEntry* m_table;
    /** \brief Total capacity of the hash table (= size of m_historyData) */
    uint m_capacity;
    //@}

    /** @name Constructors */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Default constructor */
    __HOSTDEVICE__
    ContactMemoryView()
        : m_historyData(nullptr)
        , m_table(nullptr)
        , m_capacity(0)
    {
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Constructor with parameters */
    __HOSTDEVICE__
    ContactMemoryView(ContactHistory<T>* historyData, ContactEntry* table, uint capacity)
        : m_historyData(historyData)
        , m_table(table)
        , m_capacity(capacity)
    {
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Claims a slot for a new contact, resets its history, and returns the slot index. */
    __HOSTDEVICE__
    bool claimSlot(uint slot, uint expectedState, uint2 key, uint& index)
    {
        ContactEntry& entry = m_table[slot];

#ifdef __CUDA_ARCH__
        if(atomicCAS(&entry.m_valid, expectedState, ContactEntry::ACTIVE) != expectedState)
            return false;

        m_historyData[slot] = ContactHistory<T>();
        entry.m_key         = key;
        __threadfence();
        index = slot;
        return true;
#else
        if(entry.m_valid != expectedState)
            return false;

        m_historyData[slot] = ContactHistory<T>();
        entry.m_key         = key;
        entry.m_valid       = ContactEntry::ACTIVE;
        index               = slot;
        return true;
#endif
    }
    //@}

    /** @name Hash table operations */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Finds an existing contact index in the hash table
        @param key pair of component IDs (i,j) where i < j
        @param index output parameter for the contact state index
        @return true if found, false otherwise */
    __HOSTDEVICE__
    bool find(uint2 key, uint& index) const
    {
        if(m_capacity == 0 || m_table == nullptr)
            return false;

        uint h = primeHash(key) % m_capacity;

        // Linear probing
        for(uint i = 0; i < m_capacity; ++i)
        {
            uint                idx = (h + i) % m_capacity;
            const ContactEntry& e   = m_table[idx];

            // Empty slot found - key not in table
            if(e.m_valid == ContactEntry::EMPTY)
                return false;

            // Key found — return slot position as index
            if(e.m_valid == ContactEntry::ACTIVE && e.m_key.x == key.x && e.m_key.y == key.y)
            {
#ifdef __CUDA_ARCH__
                __threadfence();
#endif
                index = idx;
                return true;
            }
        }

        // Table full, key not found
        return false;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Finds an existing contact or inserts a new one
        @param key pair of component IDs (i,j) where i < j
        @param index output parameter for the contact state index
        @return true if found or inserted, false if table is full */
    __HOSTDEVICE__
    bool findOrInsert(uint2 key, uint& index)
    {
        if(m_capacity == 0 || m_table == nullptr)
            return false;

        uint h = primeHash(key) % m_capacity;

        // Use m_capacity as sentinel: m_capacity is never a valid slot index.
        const uint NO_CLAIMABLE   = m_capacity;
        uint       firstClaimable = NO_CLAIMABLE;
        bool       retryProbe     = false;

        do
        {
            retryProbe     = false;
            firstClaimable = NO_CLAIMABLE;

            for(uint i = 0; i < m_capacity; ++i)
            {
                uint          idx = (h + i) % m_capacity;
                ContactEntry& e   = m_table[idx];

#ifdef __CUDA_ARCH__
                if(e.m_valid == ContactEntry::EMPTY)
                {
                    const uint targetSlot = firstClaimable != NO_CLAIMABLE ? firstClaimable : idx;
                    const uint expectedState = firstClaimable != NO_CLAIMABLE
                                                   ? ContactEntry::TOMBSTONE
                                                   : ContactEntry::EMPTY;
                    if(claimSlot(targetSlot, expectedState, key, index))
                        return true;

                    retryProbe = true;
                    break;
                }
                else if(e.m_valid == ContactEntry::ACTIVE && e.m_key.x == key.x
                        && e.m_key.y == key.y)
                {
                    __threadfence();
                    index = idx;
                    return true;
                }
                else if(e.m_valid == ContactEntry::TOMBSTONE && e.m_key.x == key.x
                        && e.m_key.y == key.y)
                {
                    // Contact continued — atomically promote TOMBSTONE → ACTIVE.
                    // Use CAS (not unconditional atomicExch) so we can detect the race where
                    // a concurrent thread's claimSlot() claims this same tombstone for a
                    // different key between our key-match read and the promotion write.
                    if(atomicCAS(&e.m_valid, ContactEntry::TOMBSTONE, ContactEntry::ACTIVE)
                       == ContactEntry::TOMBSTONE)
                    {
                        __threadfence();
                        index = idx;
                        return true;
                    }
                    // Another thread claimed this slot; restart the probe.
                    retryProbe = true;
                    break;
                }
                else if(e.m_valid == ContactEntry::TOMBSTONE && firstClaimable == NO_CLAIMABLE)
                {
                    firstClaimable = idx;
                }
#else
                if(e.m_valid == ContactEntry::EMPTY)
                {
                    const uint targetSlot = firstClaimable != NO_CLAIMABLE ? firstClaimable : idx;
                    const uint expectedState = firstClaimable != NO_CLAIMABLE
                                                   ? ContactEntry::TOMBSTONE
                                                   : ContactEntry::EMPTY;
                    return claimSlot(targetSlot, expectedState, key, index);
                }
                else if(e.m_valid == ContactEntry::ACTIVE && e.m_key.x == key.x
                        && e.m_key.y == key.y)
                {
                    index = idx;
                    return true;
                }
                else if(e.m_valid == ContactEntry::TOMBSTONE && e.m_key.x == key.x
                        && e.m_key.y == key.y)
                {
                    // Contact continued — promote back to active, keep history
                    e.m_valid = ContactEntry::ACTIVE;
                    index     = idx;
                    return true;
                }
                else if(e.m_valid == ContactEntry::TOMBSTONE && firstClaimable == NO_CLAIMABLE)
                {
                    firstClaimable = idx;
                }
#endif
            }

            if(!retryProbe && firstClaimable != NO_CLAIMABLE)
            {
                if(claimSlot(firstClaimable, ContactEntry::TOMBSTONE, key, index))
                    return true;
#ifdef __CUDA_ARCH__
                // On the GPU the CAS inside claimSlot can lose to a concurrent thread.
                // Rather than returning false (which propagates a null historyPtr to the
                // force kernel), restart the full probe to find another available slot.
                retryProbe = true;
#endif
            }

        } while(retryProbe);

        // Hash table is full - cannot insert
        return false;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Removes a contact from the hash table
        @param key pair of component IDs (i,j) where i < j
        @return true if removed, false if not found */
    __HOSTDEVICE__
    bool remove(uint2 key)
    {
        if(m_capacity == 0 || m_table == nullptr)
            return false;

        uint h = primeHash(key) % m_capacity;

        for(uint i = 0; i < m_capacity; ++i)
        {
            uint          idx = (h + i) % m_capacity;
            ContactEntry& e   = m_table[idx];

            // Empty slot - key not found
            if(e.m_valid == ContactEntry::EMPTY)
                return false;

            // Key found — remove it
            if(e.m_valid == ContactEntry::ACTIVE && e.m_key.x == key.x && e.m_key.y == key.y)
            {
#ifdef __CUDA_ARCH__
                // Only change the state atomically.  Do NOT reset m_historyData here: that
                // write is non-atomic and would race with another thread that claims this
                // newly-tombstoned slot via claimSlot() and is already inside computeForces()
                // reading/writing the same slot.  claimSlot() always resets the history when
                // a new pair takes the slot, so no extra reset is needed.
                atomicExch(&e.m_valid, ContactEntry::TOMBSTONE);
#else
                m_historyData[idx] = ContactHistory<T>();
                e.m_valid          = ContactEntry::TOMBSTONE;
#endif
                return true;
            }
        }

        // Key not found
        return false;
    }
    //@}
};

// =================================================================================================
/** @brief Hash table manager for contact history tracking.

    This class manages both the hash table for fast lookups and the flat array of contact history
    data. It ensures both structures are properly sized and synchronized.

    @author A.Yazdani - 2026 - Construction */
// =================================================================================================
template <typename T, MemType M = MemType::HOST>
class ContactHashTable
{
private:
    /** @name Parameters */
    //@{
    /** \brief Flat array of contact history data */
    GrainsMemBuffer<ContactHistory<T>, M> m_historyData;
    /** \brief Hash table entries buffer */
    GrainsMemBuffer<ContactEntry, M> m_table;
    /** \brief Total capacity of the hash table */
    uint m_capacity;
    /** \brief Maximum number of contacts that can be stored (used for sizing) */
    uint m_maxContacts;
    //@}

public:
    /** @name Constructors */
    //@{
    /** @brief Default constructor */
    ContactHashTable();

    /** @brief Constructor with specified capacities
        @param hashCapacity number of entries in the hash table
        @param maxContacts maximum number of contacts to store */
    ContactHashTable(uint hashCapacity, uint maxContacts);

    /** @brief Destructor */
    ~ContactHashTable();
    //@}

    /** @name Get methods */
    //@{
    /** @brief Gets the pointer to the history data */
    const ContactHistory<T>* getHistoryData() const;

    /** @brief Gets the pointer to the history data (mutable) */
    ContactHistory<T>* getHistoryData();

    /** @brief Gets the pointer to the table data */
    const ContactEntry* getTable() const;

    /** @brief Gets the pointer to the table data (mutable) */
    ContactEntry* getTable();

    /** @brief Gets the capacity of the hash table */
    uint getCapacity() const;

    /** @brief Gets the maximum number of contacts */
    uint getMaxContacts() const;

    /** @brief Counts the number of ACTIVE entries in the hash table (host-side, O(capacity)) */
    uint countActive() const;

    /** @brief Gets a complete view for passing to kernels */
    ContactMemoryView<T> getView();
    //@}

    /** @name Memory management methods */
    //@{
    /** @brief Allocates memory for both hash table and history data
        @param hashCapacity number of hash table entries to allocate
        @param maxContacts maximum number of contacts to store */
    void allocate(uint hashCapacity, uint maxContacts);

    /** @brief Frees memory used by hash table and history data */
    void deallocate();

    /** @brief Clears all entries in the hash table and resets history data */
    void clear();

    /** @brief Grows the contact table to at least newMaxContacts if it is currently smaller.
        The existing history is discarded (all contacts are treated as new after a grow).
        @param newMaxContacts the minimum required number of contact slots */
    void grow(uint newMaxContacts);

    /** @brief Performs mark-and-sweep cleanup: demotes active to stale, removes stale to empty.
        Should be called at the beginning of each timestep before contact detection */
    void markAndSweep();
    //@}
};

#endif
