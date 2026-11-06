#ifndef _CONTACTTABLE_HH_
#define _CONTACTTABLE_HH_

#include "Basic.hh"
#include "GrainsMemBuffer.hh"
#include "GrainsUtils.hh"
#include "Vector3.hh"

// =================================================================================================
/** @brief Hash table entry for contact tracking.

    Each entry stores a key (pair of component IDs), an index pointing to the ContactForce in a
    flat array, and a flag for hash table management.

    @author A.Yazdani - 2026 - Construction */
// =================================================================================================
struct ContactEntry
{
    /** @name Parameters */
    //@{
    /** \brief Contact pair key (i,j) where i < j */
    IndexPair m_key;
    /** \brief Index into the ContactForce array */
    uint m_index;
    /** \brief Validity and activity flag: 0=empty, 1=stale (not touched this step), 2=active */
    uint m_valid;
    //@}

    /** @name Constructors */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Default constructor */
    __HOSTDEVICE__
    ContactEntry()
        : m_key(makeIndexPair(0, 0))
        , m_index(0xFFFFFFFF)  // Invalid index marker
        , m_valid(0)
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

    This unified structure contains everything needed for contact history tracking: hash table
    pointers for looking up contact indices, and the flat array of contact history data. Contains
    only raw pointers and primitive types, making it safe to pass to kernels by value.

    @author A.Yazdani - 2026 - Construction */
// =================================================================================================
template <typename T>
struct ContactMemoryView
{
    /** @name Parameters */
    //@{
    /** \brief Pointer to flat array of contact history data */
    ContactHistory<T>* m_historyData;
    /** \brief Pointer to hash table entries */
    ContactEntry* m_table;
    /** \brief Total capacity of the hash table */
    uint m_capacity;
    /** \brief Pointer to the next index counter */
    uint* m_nextIndex;
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
        , m_nextIndex(nullptr)
    {
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Constructor with parameters
        @param historyData pointer to contact history data array
        @param table pointer to hash table entries
        @param capacity total capacity of the hash table
        @param nextIndex pointer to the next index counter */
    __HOSTDEVICE__
    ContactMemoryView(ContactHistory<T>* historyData,
                      ContactEntry*      table,
                      uint               capacity,
                      uint*              nextIndex)
        : m_historyData(historyData)
        , m_table(table)
        , m_capacity(capacity)
        , m_nextIndex(nextIndex)
    {
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
    /** \brief Maximum number of contacts that can be stored */
    uint m_maxContacts;
    /** \brief Pointer to the next index counter (in device/host memory) */
    uint* m_nextIndex;
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

    /** @brief Gets the current next index value (host-side read only) */
    uint getNextIndex() const;

    /** @brief Gets the pointer to the next index counter */
    uint* getNextIndexPointer();

    /** @brief Gets the pointer to the next index counter (const) */
    const uint* getNextIndexPointer() const;

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

    /** @brief Resets the index counter */
    void resetIndexCounter();

    /** @brief Performs mark-and-sweep cleanup: demotes active to stale, removes stale to empty.
        Should be called at the beginning of each timestep before contact detection */
    void markAndSweep();
    //@}
};

#endif
