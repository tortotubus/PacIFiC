#include <cuda_runtime.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <type_traits>

#include "CellsFactory.hh"
#include "GrainsCUDAError.hh"

namespace
{
    template <typename T>
    struct DeviceCellSnapshot
    {
        T              minCorner[3];
        T              maxCorner[3];
        T              linkedMinCorner[3];
        T              cellSize;
        CellDimensions dimensions;
        uint           numCells;
        uint           exists;
    };

    template <typename T, CellOrdering OrderingScheme>
    __GLOBAL__ void captureCellSnapshot(Cells<T, OrderingScheme>* const* cells,
                                        DeviceCellSnapshot<T>*           snapshot)
    {
        if(blockIdx.x != 0 || threadIdx.x != 0)
            return;

        const Cells<T, OrderingScheme>* cell = cells[0];
        if(cell == nullptr)
        {
            snapshot->exists = 0;
            return;
        }

        const Vector3<T>& minCorner       = cell->getMinCorner();
        const Vector3<T>& maxCorner       = cell->getMaxCorner();
        const Vector3<T>& linkedMinCorner = cell->getMinCornerLinkedCell();
        for(uint i = 0; i < 3; ++i)
        {
            snapshot->minCorner[i]       = minCorner[i];
            snapshot->maxCorner[i]       = maxCorner[i];
            snapshot->linkedMinCorner[i] = linkedMinCorner[i];
        }
        snapshot->cellSize   = cell->getCellSize();
        snapshot->dimensions = cell->getNumCellsPerDirection();
        snapshot->numCells   = cell->getNumCells();
        snapshot->exists     = 1;
    }

    template <typename T, CellOrdering OrderingScheme>
    __GLOBAL__ void destroyCell(Cells<T, OrderingScheme>** cells)
    {
        if(blockIdx.x == 0 && threadIdx.x == 0)
        {
            delete cells[0];
            cells[0] = nullptr;
        }
    }

    template <typename T, CellOrdering OrderingScheme>
    class HostCellOwner
    {
    public:
        ~HostCellOwner()
        {
            cells.freePointedObjects();
        }

        GrainsMemBuffer<Cells<T, OrderingScheme>*, MemType::HOST> cells;
    };

    template <typename T, CellOrdering OrderingScheme>
    class DeviceCellOwner
    {
    public:
        ~DeviceCellOwner()
        {
            if(ownsObject)
            {
                destroyCell<<<1, 1>>>(cells.getData());
                cudaDeviceSynchronize();
            }
        }

        void destroy()
        {
            destroyCell<<<1, 1>>>(cells.getData());
            cudaErrCheck(cudaDeviceSynchronize());
            ownsObject = false;
        }

        GrainsMemBuffer<Cells<T, OrderingScheme>*, MemType::DEVICE> cells;
        bool                                                        ownsObject = false;
    };

    template <typename T, CellOrdering OrderingScheme>
    struct CellFactoryCase
    {
        using Real                             = T;
        static constexpr CellOrdering ordering = OrderingScheme;
    };

    using CellFactoryCases = ::testing::Types<CellFactoryCase<float, CellOrdering::LINEAR>,
                                              CellFactoryCase<double, CellOrdering::LINEAR>,
                                              CellFactoryCase<float, CellOrdering::MORTON>,
                                              CellFactoryCase<double, CellOrdering::MORTON>>;

    template <typename Case>
    class CellsFactoryDeviceTest : public ::testing::Test
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

    TYPED_TEST_SUITE(CellsFactoryDeviceTest, CellFactoryCases);

    TYPED_TEST(CellsFactoryDeviceTest, DeviceConstructionMatchesHost)
    {
        using T                               = typename TypeParam::Real;
        constexpr CellOrdering OrderingScheme = TypeParam::ordering;
        const Vector3<T>       minCorner(T(-1), T(-2), T(-3));
        const Vector3<T>       maxCorner(T(2), T(2), T(2));
        const T                requestedCellSize = T(0.75);

        HostCellOwner<T, OrderingScheme> host;
        const uint                       hostCount
            = CellsFactory<T, OrderingScheme>::template create<MemType::HOST>(minCorner,
                                                                              maxCorner,
                                                                              requestedCellSize,
                                                                              host.cells);

        DeviceCellOwner<T, OrderingScheme> device;
        const uint                         deviceCount
            = CellsFactory<T, OrderingScheme>::template create<MemType::DEVICE>(minCorner,
                                                                                maxCorner,
                                                                                requestedCellSize,
                                                                                device.cells);
        device.ownsObject = true;

        GrainsMemBuffer<DeviceCellSnapshot<T>, MemType::DEVICE> deviceSnapshot(1);
        GrainsMemBuffer<DeviceCellSnapshot<T>, MemType::HOST>   hostSnapshot;
        captureCellSnapshot<<<1, 1>>>(device.cells.getData(), deviceSnapshot.getData());
        ASSERT_EQ(cudaDeviceSynchronize(), cudaSuccess);
        hostSnapshot.copyFrom(deviceSnapshot);
        ASSERT_EQ(cudaDeviceSynchronize(), cudaSuccess);

        const Cells<T, OrderingScheme>* expected           = host.cells[0];
        const auto&                     actual             = hostSnapshot[0];
        const CellDimensions            expectedDimensions = expected->getNumCellsPerDirection();
        const T                         tolerance = std::is_same_v<T, float> ? T(1e-6) : T(1e-13);

        EXPECT_EQ(actual.exists, 1u);
        EXPECT_EQ(deviceCount, hostCount);
        EXPECT_EQ(actual.numCells, hostCount);
        EXPECT_EQ(actual.dimensions.x, expectedDimensions.x);
        EXPECT_EQ(actual.dimensions.y, expectedDimensions.y);
        EXPECT_EQ(actual.dimensions.z, expectedDimensions.z);
        EXPECT_EQ(actual.dimensions.w, expectedDimensions.w);
        EXPECT_NEAR(actual.cellSize, expected->getCellSize(), tolerance);
        for(uint i = 0; i < 3; ++i)
        {
            EXPECT_NEAR(actual.minCorner[i], expected->getMinCorner()[i], tolerance);
            EXPECT_NEAR(actual.maxCorner[i], expected->getMaxCorner()[i], tolerance);
            EXPECT_NEAR(actual.linkedMinCorner[i],
                        expected->getMinCornerLinkedCell()[i],
                        tolerance);
        }

        device.destroy();
        GrainsMemBuffer<Cells<T, OrderingScheme>*, MemType::HOST> releasedPointers;
        releasedPointers.copyFrom(device.cells);
        ASSERT_EQ(cudaDeviceSynchronize(), cudaSuccess);
        EXPECT_EQ(releasedPointers[0], nullptr);
    }
}  // namespace
