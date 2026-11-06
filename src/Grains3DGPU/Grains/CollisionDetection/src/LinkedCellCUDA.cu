#include "GrainsPlatform.hh"
#include "LinkedCellCUDA.hh"

#include <cub/device/device_radix_sort.cuh>
#include <cub/device/device_reduce.cuh>
#include <cub/device/device_scan.cuh>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>

#include "GrainsCUDAUtils.hh"
#include "GrainsParameters.hh"
#include "LinkedCell_Kernels.hh"
#include "VectorMath.hh"

namespace grains::cuda
{
    namespace detail
    {
        std::size_t exclusiveSumStorageBytes(const unsigned int* input,
                                             unsigned int*       output,
                                             unsigned int        count)
        {
            if(count == 0)
                return 0;
            std::size_t bytes = 0;
            cudaErrCheck(cub::DeviceScan::ExclusiveSum(nullptr, bytes, input, output, count));
            return bytes;
        }

        void exclusiveSum(const unsigned int* input,
                          unsigned int*       output,
                          unsigned int        count,
                          void*               scratch,
                          std::size_t         scratchBytes)
        {
            if(count == 0)
                return;
            cudaErrCheck(
                cub::DeviceScan::ExclusiveSum(scratch, scratchBytes, input, output, count));
        }

        void exclusiveSumAsync(const unsigned int* input,
                               unsigned int*       output,
                               unsigned int        count,
                               void*               scratch,
                               std::size_t         scratchBytes,
                               void*               stream)
        {
            if(count == 0)
                return;
            cudaErrCheck(cub::DeviceScan::ExclusiveSum(scratch,
                                                       scratchBytes,
                                                       input,
                                                       output,
                                                       count,
                                                       static_cast<cudaStream_t>(stream)));
        }

        std::size_t sortCellKeysStorageBytes(const std::uint64_t* input,
                                             std::uint64_t*       output,
                                             unsigned int         count)
        {
            if(count == 0)
                return 0;
            std::size_t bytes = 0;
            cudaErrCheck(cub::DeviceRadixSort::SortKeys(nullptr, bytes, input, output, count));
            return bytes;
        }

        void sortCellKeys(const std::uint64_t* input,
                          std::uint64_t*       output,
                          unsigned int         count,
                          void*                scratch,
                          std::size_t          scratchBytes)
        {
            if(count == 0)
                return;
            cudaErrCheck(
                cub::DeviceRadixSort::SortKeys(scratch, scratchBytes, input, output, count));
        }
    }

    struct SortBasedState
    {
        GrainsMemBuffer<std::uint64_t, MemType::DEVICE> cellParticleIDs;
        GrainsMemBuffer<unsigned int, MemType::DEVICE>  cellPrefixSums;
        void*                                           sortTemporaryStorage      = nullptr;
        std::size_t                                     sortTemporaryStorageBytes = 0;
        unsigned int*                                   deviceNumCells            = nullptr;
        cudaEvent_t                                     resizeComplete            = nullptr;
        cudaStream_t                                    neighborStream            = nullptr;
        cudaStream_t                                    oldPositionStream         = nullptr;
        cudaStream_t                                    particleStream            = nullptr;

        SortBasedState(unsigned int numParticles, unsigned int numCells)
            : cellParticleIDs(numParticles)
            , cellPrefixSums(numCells)
        {
        }

        ~SortBasedState()
        {
            if(sortTemporaryStorage != nullptr)
                cudaFree(sortTemporaryStorage);
            if(deviceNumCells != nullptr)
                cudaFree(deviceNumCells);
            if(resizeComplete != nullptr)
                cudaEventDestroy(resizeComplete);
            if(neighborStream != nullptr)
                cudaStreamDestroy(neighborStream);
            if(oldPositionStream != nullptr)
                cudaStreamDestroy(oldPositionStream);
            if(particleStream != nullptr)
                cudaStreamDestroy(particleStream);
        }
    };

    void SortBasedStateDeleter::operator()(SortBasedState* state) const
    {
        delete state;
    }

    SortBasedStatePtr createSortBasedState(unsigned int numParticles, unsigned int numCells)
    {
        SortBasedStatePtr state(new SortBasedState(numParticles, numCells));
        state->sortTemporaryStorageBytes
            = detail::sortCellKeysStorageBytes(state->cellParticleIDs.getData(),
                                               state->cellParticleIDs.getData(),
                                               numParticles);
        if(state->sortTemporaryStorageBytes != 0)
            cudaErrCheck(
                cudaMalloc(&state->sortTemporaryStorage, state->sortTemporaryStorageBytes));
        cudaErrCheck(cudaMalloc(&state->deviceNumCells, sizeof(unsigned int)));
        cudaErrCheck(cudaEventCreate(&state->resizeComplete));
        cudaErrCheck(cudaStreamCreate(&state->neighborStream));
        cudaErrCheck(cudaStreamCreate(&state->oldPositionStream));
        cudaErrCheck(cudaStreamCreate(&state->particleStream));
        return state;
    }

    const std::uint64_t* sortBasedCellParticleIDs(const SortBasedState& state)
    {
        return state.cellParticleIDs.getData();
    }

    const unsigned int* sortBasedCellPrefixSums(const SortBasedState& state)
    {
        return state.cellPrefixSums.getData();
    }

    template <typename T>
    unsigned int prepareSortBasedUpdate(SortBasedState& state, SortBasedUpdateView<T> view)
    {
        unsigned int numBlocks;
        unsigned int numThreads;
        computeOptimalThreadsAndBlocks(view.numParticles,
                                       GrainsParameters<T>::m_GPU,
                                       numBlocks,
                                       numThreads);

        unsigned int numCells = view.numCells;
        if(view.adaptive)
        {
            resizeCells_Device<<<1, 1>>>(view.cells, view.cellSize, state.deviceNumCells);
            cudaErrCheck(cudaGetLastError());
            cudaErrCheck(cudaMemcpyAsync(&numCells,
                                         state.deviceNumCells,
                                         sizeof(unsigned int),
                                         cudaMemcpyDeviceToHost,
                                         nullptr));
            cudaErrCheck(cudaStreamSynchronize(nullptr));
            cudaErrCheck(cudaEventRecord(state.resizeComplete, nullptr));

            view.neighborCells->reserve(numCells * 27, state.neighborStream);
            state.cellPrefixSums.reserve(numCells, state.particleStream);

            cudaErrCheck(
                cudaMemcpyAsync(view.oldPositions,
                                view.positions,
                                (view.numObstacles + view.numParticles) * sizeof(Vector3<T>),
                                cudaMemcpyDeviceToDevice,
                                state.oldPositionStream));

            unsigned int numCellBlocks;
            unsigned int numCellThreads;
            computeOptimalThreadsAndBlocks(numCells,
                                           GrainsParameters<T>::m_GPU,
                                           numCellBlocks,
                                           numCellThreads);
            cudaErrCheck(cudaStreamWaitEvent(state.neighborStream, state.resizeComplete, 0));
            view.neighborCells->fill(UINT_MAX, state.neighborStream);
            generateNeighborCells_Device<<<numCellBlocks,
                                           numCellThreads,
                                           0,
                                           state.neighborStream>>>(view.cells,
                                                                   numCells,
                                                                   view.neighborCells->getData());
            cudaErrCheck(cudaGetLastError());
        }
        else
        {
            // Preserve the ordering guarantee previously provided by the null-stream sync.
            cudaErrCheck(cudaStreamSynchronize(nullptr));
        }

        state.cellPrefixSums.fill(UINT_MAX, state.particleStream);
        computeCellParticleIDs_Device<<<numBlocks, numThreads, 0, state.particleStream>>>(
            view.cells,
            view.positions + view.numObstacles,
            view.numParticles,
            view.numObstacles,
            state.cellParticleIDs.getData());
        cudaErrCheck(cudaGetLastError());

        if(view.adaptive)
        {
            cudaErrCheck(cudaStreamSynchronize(state.neighborStream));
            cudaErrCheck(cudaStreamSynchronize(state.oldPositionStream));
        }
        cudaErrCheck(cudaStreamSynchronize(state.particleStream));
        return numCells;
    }

    template <typename T>
    void finalizeSortBasedUpdate(SortBasedState& state, unsigned int numParticles)
    {
        detail::sortCellKeys(state.cellParticleIDs.getData(),
                             state.cellParticleIDs.getData(),
                             numParticles,
                             state.sortTemporaryStorage,
                             state.sortTemporaryStorageBytes);

        unsigned int numBlocks;
        unsigned int numThreads;
        computeOptimalThreadsAndBlocks(numParticles,
                                       GrainsParameters<T>::m_GPU,
                                       numBlocks,
                                       numThreads);
        const unsigned int sharedMemoryBytes = sizeof(unsigned int) * (numThreads + 1);
        computeCellStart_Kernel<<<numBlocks, numThreads, sharedMemoryBytes>>>(
            state.cellParticleIDs.getData(),
            numParticles,
            state.cellPrefixSums.getData());
        cudaErrCheck(cudaGetLastError());
        cudaErrCheck(cudaStreamSynchronize(nullptr));
    }

    template <typename T>
    struct LinkedCellReductionState
    {
        void*        temporaryStorage      = nullptr;
        std::size_t  temporaryStorageBytes = 0;
        cudaStream_t stream                = nullptr;

        ~LinkedCellReductionState()
        {
            if(temporaryStorage != nullptr)
                cudaFree(temporaryStorage);
            if(stream != nullptr)
                cudaStreamDestroy(stream);
        }
    };

    namespace
    {
        template <typename T>
        struct DisplacementSquared
        {
            const Vector3<T>* oldPositions;
            const Vector3<T>* positions;

            __device__ T operator()(std::size_t i) const
            {
                return norm2(positions[i] - oldPositions[i]);
            }
        };

        template <typename T>
        cudaError_t reduceMaxDisplacementSquared(void*             temporaryStorage,
                                                 std::size_t&      temporaryStorageBytes,
                                                 const Vector3<T>* oldPositions,
                                                 const Vector3<T>* positions,
                                                 std::size_t       count,
                                                 T*                output,
                                                 cudaStream_t      stream)
        {
            auto indices = thrust::make_counting_iterator<std::size_t>(0);
            auto input
                = thrust::make_transform_iterator(indices,
                                                  DisplacementSquared<T>{oldPositions, positions});
            return cub::DeviceReduce::Max(temporaryStorage,
                                          temporaryStorageBytes,
                                          input,
                                          output,
                                          count,
                                          stream);
        }
    }

    void LinkedCellReductionStateDeleter::operator()(LinkedCellReductionState<float>* state) const
    {
        delete state;
    }

    void LinkedCellReductionStateDeleter::operator()(LinkedCellReductionState<double>* state) const
    {
        delete state;
    }

    template <typename T>
    LinkedCellReductionStatePtr<T> createLinkedCellReductionState(std::size_t       capacity,
                                                                  const Vector3<T>* oldPositions,
                                                                  const Vector3<T>* positions,
                                                                  T*                output)
    {
        LinkedCellReductionStatePtr<T> state(new LinkedCellReductionState<T>);
        cudaErrCheck(cudaStreamCreate(&state->stream));
        cudaErrCheck(reduceMaxDisplacementSquared(state->temporaryStorage,
                                                  state->temporaryStorageBytes,
                                                  oldPositions,
                                                  positions,
                                                  capacity,
                                                  output,
                                                  state->stream));
        if(state->temporaryStorageBytes != 0)
            cudaErrCheck(cudaMalloc(&state->temporaryStorage, state->temporaryStorageBytes));
        return state;
    }

    template <typename T>
    void computeMaxDisplacementSquared(LinkedCellReductionState<T>& state,
                                       const Vector3<T>*            oldPositions,
                                       const Vector3<T>*            positions,
                                       std::size_t                  count,
                                       T*                           output)
    {
        cudaErrCheck(reduceMaxDisplacementSquared(state.temporaryStorage,
                                                  state.temporaryStorageBytes,
                                                  oldPositions,
                                                  positions,
                                                  count,
                                                  output,
                                                  state.stream));
        cudaErrCheck(cudaStreamSynchronize(state.stream));
    }

    template <typename T>
    void launchMaxDisplacementSquared(LinkedCellReductionState<T>& state,
                                      const Vector3<T>*            oldPositions,
                                      const Vector3<T>*            positions,
                                      std::size_t                  count,
                                      T*                           output)
    {
        cudaErrCheck(reduceMaxDisplacementSquared(state.temporaryStorage,
                                                  state.temporaryStorageBytes,
                                                  oldPositions,
                                                  positions,
                                                  count,
                                                  output,
                                                  state.stream));
    }

    template <typename T>
    void synchronizeLinkedCellReduction(LinkedCellReductionState<T>& state)
    {
        cudaErrCheck(cudaStreamSynchronize(state.stream));
    }

    template <typename T>
    void generateNeighborCells(const Cells<T>* const* cells,
                               unsigned int           numCells,
                               unsigned int*          neighborCells)
    {
        unsigned int numBlocks;
        unsigned int numThreads;
        computeOptimalThreadsAndBlocks(numCells, GrainsParameters<T>::m_GPU, numBlocks, numThreads);
        generateNeighborCells_Device<<<numBlocks, numThreads>>>(cells, numCells, neighborCells);
        cudaErrCheck(cudaGetLastError());
        cudaErrCheck(cudaDeviceSynchronize());
    }

    template <typename T>
    void linkObstacles(const RigidBody<T>* const* rigidBodies,
                       const Vector3<T>*          positions,
                       const Quaternion<T>*       quaternions,
                       const Cells<T>* const*     cells,
                       unsigned int               numObstacles,
                       unsigned int               maxCellsPerObstacle,
                       IndexPair*                 obstacleIDs,
                       unsigned int*              obstacleCellIDs)
    {
        linkObstacles_Device<T><<<numObstacles, 1>>>(rigidBodies,
                                                     positions,
                                                     quaternions,
                                                     cells,
                                                     numObstacles,
                                                     maxCellsPerObstacle,
                                                     obstacleIDs,
                                                     obstacleCellIDs);
        cudaErrCheck(cudaGetLastError());
        cudaErrCheck(cudaDeviceSynchronize());
    }

    template LinkedCellReductionStatePtr<float> createLinkedCellReductionState<float>(
        std::size_t, const Vector3<float>*, const Vector3<float>*, float*);
    template LinkedCellReductionStatePtr<double> createLinkedCellReductionState<double>(
        std::size_t, const Vector3<double>*, const Vector3<double>*, double*);
    template void         computeMaxDisplacementSquared<float>(LinkedCellReductionState<float>&,
                                                       const Vector3<float>*,
                                                       const Vector3<float>*,
                                                       std::size_t,
                                                       float*);
    template void         computeMaxDisplacementSquared<double>(LinkedCellReductionState<double>&,
                                                        const Vector3<double>*,
                                                        const Vector3<double>*,
                                                        std::size_t,
                                                        double*);
    template void         launchMaxDisplacementSquared<float>(LinkedCellReductionState<float>&,
                                                      const Vector3<float>*,
                                                      const Vector3<float>*,
                                                      std::size_t,
                                                      float*);
    template void         launchMaxDisplacementSquared<double>(LinkedCellReductionState<double>&,
                                                       const Vector3<double>*,
                                                       const Vector3<double>*,
                                                       std::size_t,
                                                       double*);
    template void         synchronizeLinkedCellReduction<float>(LinkedCellReductionState<float>&);
    template void         synchronizeLinkedCellReduction<double>(LinkedCellReductionState<double>&);
    template unsigned int prepareSortBasedUpdate<float>(SortBasedState&,
                                                        SortBasedUpdateView<float>);
    template unsigned int prepareSortBasedUpdate<double>(SortBasedState&,
                                                         SortBasedUpdateView<double>);
    template void         finalizeSortBasedUpdate<float>(SortBasedState&, unsigned int);
    template void         finalizeSortBasedUpdate<double>(SortBasedState&, unsigned int);
    template void
        generateNeighborCells<float>(const Cells<float>* const*, unsigned int, unsigned int*);
    template void
        generateNeighborCells<double>(const Cells<double>* const*, unsigned int, unsigned int*);
    template void linkObstacles<float>(const RigidBody<float>* const*,
                                       const Vector3<float>*,
                                       const Quaternion<float>*,
                                       const Cells<float>* const*,
                                       unsigned int,
                                       unsigned int,
                                       IndexPair*,
                                       unsigned int*);
    template void linkObstacles<double>(const RigidBody<double>* const*,
                                        const Vector3<double>*,
                                        const Quaternion<double>*,
                                        const Cells<double>* const*,
                                        unsigned int,
                                        unsigned int,
                                        IndexPair*,
                                        unsigned int*);
}
