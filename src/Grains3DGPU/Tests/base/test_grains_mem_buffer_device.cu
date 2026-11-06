#include "GrainsMemBuffer.hh"

#include <cuda_runtime.h>
#include <gtest/gtest.h>

#include <limits>

namespace
{
    class GrainsMemBufferDeviceTest : public ::testing::Test
    {
    protected:
        void SetUp() override
        {
            int               deviceCount = 0;
            const cudaError_t error       = cudaGetDeviceCount(&deviceCount);
            if(error != cudaSuccess || deviceCount == 0)
                GTEST_SKIP() << "No CUDA device available";
        }
    };

    TEST_F(GrainsMemBufferDeviceTest, CopiesFillsSequencesAndReallocatesOnAStream)
    {
        constexpr size_t count = 8;
        cudaStream_t     stream;
        ASSERT_EQ(cudaStreamCreate(&stream), cudaSuccess);

        GrainsMemBuffer<unsigned int, MemType::DEVICE> device;
        device.initialize(count);
        device.fill(7u, stream);

        GrainsMemBuffer<unsigned int> host(count);
        device.copyTo(host, stream);
        ASSERT_EQ(cudaStreamSynchronize(stream), cudaSuccess);
        for(const unsigned int value : host)
            EXPECT_EQ(value, 7u);

        device.fill(std::numeric_limits<unsigned int>::max(), stream);
        device.copyTo(host, stream);
        ASSERT_EQ(cudaStreamSynchronize(stream), cudaSuccess);
        for(const unsigned int value : host)
            EXPECT_EQ(value, std::numeric_limits<unsigned int>::max());

        device.fill(stream);
        device.copyTo(host, stream);
        ASSERT_EQ(cudaStreamSynchronize(stream), cudaSuccess);
        for(const unsigned int value : host)
            EXPECT_EQ(value, 0u);

        device.sequence(5u, stream);
        device.reserve(16, stream);
        EXPECT_EQ(device.getSize(), count);
        EXPECT_EQ(device.getCapacity(), 16u);

        device.copyTo(host, stream);
        ASSERT_EQ(cudaStreamSynchronize(stream), cudaSuccess);
        for(size_t i = 0; i < count; ++i)
            EXPECT_EQ(host[i], static_cast<unsigned int>(i + 5));

        device.shrink_to_fit(stream);
        EXPECT_EQ(device.getCapacity(), count);
        device.copyTo(host, stream);
        ASSERT_EQ(cudaStreamSynchronize(stream), cudaSuccess);
        for(size_t i = 0; i < count; ++i)
            EXPECT_EQ(host[i], static_cast<unsigned int>(i + 5));

        ASSERT_EQ(cudaStreamDestroy(stream), cudaSuccess);
    }

    TEST_F(GrainsMemBufferDeviceTest, SupportsFloatingPointScalarOperations)
    {
        constexpr size_t                         count = 6;
        GrainsMemBuffer<double, MemType::DEVICE> device(count, 2.5);
        GrainsMemBuffer<double>                  host(count);

        device.copyTo(host);
        ASSERT_EQ(cudaDeviceSynchronize(), cudaSuccess);
        for(const double value : host)
            EXPECT_DOUBLE_EQ(value, 2.5);

        device.sequence(1.25);
        device.copyTo(host);
        ASSERT_EQ(cudaDeviceSynchronize(), cudaSuccess);
        for(size_t i = 0; i < count; ++i)
            EXPECT_DOUBLE_EQ(host[i], 1.25 + static_cast<double>(i));
    }

    TEST_F(GrainsMemBufferDeviceTest, AllocatesManagedPinnedAndMappedMemory)
    {
        GrainsMemBuffer<int, MemType::MANAGED> managed;
        managed.initialize(3);
        managed[0] = 1;
        managed[1] = 2;
        managed[2] = 3;
        EXPECT_EQ(managed.getDeviceData(), managed.getData());

        GrainsMemBuffer<int, MemType::PINNED> pinned;
        pinned.initialize(2);
        pinned[0] = 4;
        pinned[1] = 5;
        EXPECT_EQ(pinned.getDeviceData(), pinned.getData());

        GrainsMemBuffer<int, MemType::MAPPED> mapped;
        mapped.initialize(2);
        mapped[0] = 6;
        mapped[1] = 7;
        EXPECT_NE(mapped.getDeviceData(), nullptr);
    }

    TEST_F(GrainsMemBufferDeviceTest, FreesObjectsOwnedByDevicePointerBuffer)
    {
        int* first  = nullptr;
        int* second = nullptr;
        ASSERT_EQ(cudaMalloc(&first, sizeof(int)), cudaSuccess);
        ASSERT_EQ(cudaMalloc(&second, sizeof(int)), cudaSuccess);

        GrainsMemBuffer<int*> hostPointers(2);
        hostPointers[0] = first;
        hostPointers[1] = second;

        GrainsMemBuffer<int*, MemType::DEVICE> devicePointers(2);
        hostPointers.copyTo(devicePointers);
        ASSERT_EQ(cudaDeviceSynchronize(), cudaSuccess);

        devicePointers.freePointedObjects();
        devicePointers.copyTo(hostPointers);
        ASSERT_EQ(cudaDeviceSynchronize(), cudaSuccess);
        EXPECT_EQ(hostPointers[0], nullptr);
        EXPECT_EQ(hostPointers[1], nullptr);
    }
}  // namespace
