#include "GrainsMemBuffer.hh"

#include <gtest/gtest.h>

#include <utility>

TEST(GrainsMemBufferHostTest, GrowsAndShrinksWithoutChangingLogicalContents)
{
    GrainsMemBuffer<int> buffer(4);
    buffer.sequence(3);

    buffer.reserve(12);
    EXPECT_EQ(buffer.getSize(), 4u);
    EXPECT_EQ(buffer.getCapacity(), 12u);
    for(size_t i = 0; i < buffer.getSize(); ++i)
        EXPECT_EQ(buffer[i], static_cast<int>(i + 3));

    buffer.shrink_to_fit();
    EXPECT_EQ(buffer.getSize(), 4u);
    EXPECT_EQ(buffer.getCapacity(), 4u);
    for(size_t i = 0; i < buffer.getSize(); ++i)
        EXPECT_EQ(buffer[i], static_cast<int>(i + 3));
}

TEST(GrainsMemBufferHostTest, SupportsCopyMoveFillAndBulkInsertion)
{
    GrainsMemBuffer<unsigned int> source(3, 7u);
    const unsigned int            appended[] = {11u, 13u};
    source.push_bulk(appended, 2);

    GrainsMemBuffer<unsigned int> copy(source);
    ASSERT_EQ(copy.getSize(), 5u);
    EXPECT_EQ(copy[0], 7u);
    EXPECT_EQ(copy[2], 7u);
    EXPECT_EQ(copy[3], 11u);
    EXPECT_EQ(copy[4], 13u);

    GrainsMemBuffer<unsigned int> moved(std::move(copy));
    EXPECT_EQ(copy.getData(), nullptr);
    EXPECT_EQ(copy.getSize(), 0u);
    ASSERT_EQ(moved.getSize(), 5u);

    moved.fill(5u);
    for(const unsigned int value : moved)
        EXPECT_EQ(value, 5u);
}

TEST(GrainsMemBufferHostTest, DeletesObjectsOwnedByPointerBuffer)
{
    GrainsMemBuffer<int*> pointers(2);
    pointers[0] = new int(17);
    pointers[1] = new int(23);

    pointers.freePointedObjects();

    EXPECT_EQ(pointers[0], nullptr);
    EXPECT_EQ(pointers[1], nullptr);
}
