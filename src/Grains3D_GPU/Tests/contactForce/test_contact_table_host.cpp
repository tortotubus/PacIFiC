#include <gtest/gtest.h>

#include <array>

#include "ContactTable.hh"
#include "GrainsUtils.hh"

// Host-side tests for ContactHashTable.
// find / findOrInsert / remove live on ContactMemoryView<T> (the lightweight proxy
// returned by getView()), so each test that needs lookup operations uses a view.
// The raw pointers in the view remain valid for the lifetime of the table; a single
// view obtained in SetUp() stays usable after clear().

class ContactTableHostTest : public ::testing::Test
{
protected:
    ContactHashTable<float, MemType::HOST> table;
    ContactMemoryView<float>               view;

    void SetUp() override
    {
        table.allocate(100, 100);
        view = table.getView();
    }

    void TearDown() override
    {
        table.deallocate();
    }
};

TEST_F(ContactTableHostTest, AllocationAndDeallocation)
{
    ContactHashTable<float, MemType::HOST> t;
    EXPECT_EQ(t.getCapacity(), 0);
    EXPECT_EQ(t.countActive(), 0u);
    t.allocate(100, 100);
    EXPECT_EQ(t.getCapacity(), 100);
    EXPECT_EQ(t.countActive(), 0u);
    t.deallocate();
    EXPECT_EQ(t.getCapacity(), 0);
    EXPECT_EQ(t.countActive(), 0u);
}

TEST_F(ContactTableHostTest, ConstructorWithCapacity)
{
    ContactHashTable<float, MemType::HOST> t(50, 50);
    EXPECT_EQ(t.getCapacity(), 50);
    EXPECT_EQ(t.countActive(), 0u);
}

TEST_F(ContactTableHostTest, BasicInsertion)
{
    uint2 pair1 = make_uint2(5, 10);
    uint  index1;
    bool  success = view.findOrInsert(pair1, index1);
    EXPECT_TRUE(success);
    EXPECT_LT(index1, 100u);
    EXPECT_EQ(table.countActive(), 1u);
}

TEST_F(ContactTableHostTest, FindExisting)
{
    uint2 pair1 = make_uint2(5, 10);
    uint  index1, index2;
    view.findOrInsert(pair1, index1);
    bool found = view.find(pair1, index2);
    EXPECT_TRUE(found);
    EXPECT_EQ(index1, index2);
}

TEST_F(ContactTableHostTest, FindNonExistent)
{
    uint2 pair1 = make_uint2(5, 10);
    uint2 pair2 = make_uint2(7, 12);
    uint  index;
    view.findOrInsert(pair1, index);
    bool found = view.find(pair2, index);
    EXPECT_FALSE(found);
}

TEST_F(ContactTableHostTest, MultipleInsertions)
{
    const int numPairs = 10;
    uint      indices[numPairs];
    for(int i = 0; i < numPairs; ++i)
    {
        uint2 pair    = make_uint2(i, i + 10);
        bool  success = view.findOrInsert(pair, indices[i]);
        EXPECT_TRUE(success);
        EXPECT_LT(indices[i], 100u);
    }
    EXPECT_EQ(table.countActive(), (uint)numPairs);
}

TEST_F(ContactTableHostTest, DuplicateInsertion)
{
    uint2 pair = make_uint2(5, 10);
    uint  index1, index2;
    view.findOrInsert(pair, index1);
    view.findOrInsert(pair, index2);
    EXPECT_EQ(index1, index2);
    EXPECT_EQ(table.countActive(), 1u);
}

TEST_F(ContactTableHostTest, Removal)
{
    uint2 pair = make_uint2(5, 10);
    uint  index;
    view.findOrInsert(pair, index);
    bool removed = view.remove(pair);
    EXPECT_TRUE(removed);
    bool found = view.find(pair, index);
    EXPECT_FALSE(found);
}

TEST_F(ContactTableHostTest, RemoveNonExistent)
{
    uint2 pair    = make_uint2(5, 10);
    bool  removed = view.remove(pair);
    EXPECT_FALSE(removed);
}

TEST_F(ContactTableHostTest, ClearOperation)
{
    for(int i = 0; i < 5; ++i)
    {
        uint2 pair = make_uint2(i, i + 10);
        uint  index;
        view.findOrInsert(pair, index);
    }
    EXPECT_EQ(table.countActive(), 5u);
    table.clear();
    EXPECT_EQ(table.countActive(), 0u);
    uint index;
    for(int i = 0; i < 5; ++i)
    {
        uint2 pair  = make_uint2(i, i + 10);
        bool  found = view.find(pair, index);
        EXPECT_FALSE(found);
    }
}

TEST_F(ContactTableHostTest, FullTable)
{
    const int                              capacity = 10;
    ContactHashTable<float, MemType::HOST> t(capacity, capacity);
    ContactMemoryView<float>               localView = t.getView();
    for(int i = 0; i < capacity; ++i)
    {
        uint2 pair = make_uint2(i, i + 100);
        uint  index;
        bool  success = localView.findOrInsert(pair, index);
        EXPECT_TRUE(success);
    }
}

TEST_F(ContactTableHostTest, HashCollisions)
{
    uint indices[5];
    for(int i = 0; i < 5; ++i)
    {
        uint2 pair    = make_uint2(i * 1000, i * 1000 + 1);
        bool  success = view.findOrInsert(pair, indices[i]);
        EXPECT_TRUE(success);
    }
    for(int i = 0; i < 5; ++i)
    {
        for(int j = i + 1; j < 5; ++j)
            EXPECT_NE(indices[i], indices[j]);
    }
}

TEST_F(ContactTableHostTest, PairOrdering)
{
    ContactHashTable<float, MemType::HOST> t(100, 100);
    ContactMemoryView<float>               localView = t.getView();
    uint2                                  pair1     = make_uint2(5, 10);
    uint2                                  pair2     = make_uint2(10, 5);
    uint                                   index1, index2;
    localView.findOrInsert(pair1, index1);
    localView.findOrInsert(pair2, index2);
    EXPECT_NE(index1, index2);
}

TEST_F(ContactTableHostTest, CleanupReusesFreedHistorySlots)
{
    ContactHashTable<float, MemType::HOST> t(16, 4);
    ContactMemoryView<float>               localView = t.getView();
    std::array<uint, 4>                    recycled{};

    for(uint i = 0; i < 4; ++i)
    {
        uint index;
        bool success = localView.findOrInsert(make_uint2(i, i + 10), index);
        EXPECT_TRUE(success);
        EXPECT_LT(index, 16u);
    }

    EXPECT_EQ(t.countActive(), 4u);

    t.markAndSweep();
    t.markAndSweep();

    for(uint i = 0; i < 4; ++i)
    {
        bool success = localView.findOrInsert(make_uint2(i + 100, i + 200), recycled[i]);
        EXPECT_TRUE(success);
        EXPECT_LT(recycled[i], 16u);
    }

    EXPECT_EQ(t.countActive(), 4u);
}

TEST_F(ContactTableHostTest, CleanupKeepsProbeChainsSearchable)
{
    const uint capacity = 8u;
    uint2      pairA    = make_uint2(0u, 0u);
    uint2      pairB    = make_uint2(0u, 0u);
    bool       found    = false;

    for(uint i = 0; i < 32u && !found; ++i)
    {
        for(uint j = i + 1; j < 32u && !found; ++j)
        {
            const uint2 candidateA = make_uint2(i, j);
            for(uint k = 0; k < 32u && !found; ++k)
            {
                for(uint l = k + 1; l < 32u && !found; ++l)
                {
                    const uint2 candidateB = make_uint2(k, l);
                    if((candidateA.x != candidateB.x || candidateA.y != candidateB.y)
                       && primeHash(candidateA) % capacity == primeHash(candidateB) % capacity)
                    {
                        pairA = candidateA;
                        pairB = candidateB;
                        found = true;
                    }
                }
            }
        }
    }

    ASSERT_TRUE(found);

    ContactHashTable<float, MemType::HOST> t(capacity, capacity);
    ContactMemoryView<float>               localView = t.getView();

    uint indexA, indexB, reactivatedIndex;
    ASSERT_TRUE(localView.findOrInsert(pairA, indexA));
    ASSERT_TRUE(localView.findOrInsert(pairB, indexB));

    t.markAndSweep();
    ASSERT_TRUE(localView.findOrInsert(pairB, reactivatedIndex));
    EXPECT_EQ(reactivatedIndex, indexB);

    t.markAndSweep();

    uint foundIndex;
    EXPECT_TRUE(localView.find(pairB, foundIndex));
    EXPECT_EQ(foundIndex, indexB);
    EXPECT_EQ(t.countActive(), 2u);
}
