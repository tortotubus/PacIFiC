#include "NeighborListLinkedCellCUDA.hh"

#include "GrainsCUDAUtils.hh"
#include "GrainsParameters.hh"
#include "GrainsUtils.hh"
#include "LinkedCell.hh"
#include "LinkedCellCUDA.hh"
#include "LinkedCell_AtomicFixed.cuh"
#include "NeighborList_Kernels.cuh"

namespace grains::cuda
{
    struct NeighborListLinkedCellState
    {
        GrainsMemBuffer<unsigned int, MemType::DEVICE> numNeighbors;
        GrainsMemBuffer<unsigned int, MemType::DEVICE> numNeighborsPrefixSums;
        void*                                          scanTemporaryStorage      = nullptr;
        std::size_t                                    scanTemporaryStorageBytes = 0;
        unsigned int*                                  obstacleParticlePairCount = nullptr;
        cudaStream_t                                   obstacleStream            = nullptr;
        cudaStream_t                                   particleStream            = nullptr;

        explicit NeighborListLinkedCellState(unsigned int numParticles)
            : numNeighbors(numParticles)
            , numNeighborsPrefixSums(numParticles + 1)
        {
            numNeighbors.fill();
            numNeighborsPrefixSums.fill();

            scanTemporaryStorageBytes
                = detail::exclusiveSumStorageBytes(numNeighbors.getData(),
                                                   numNeighborsPrefixSums.getData(),
                                                   numParticles);
            if(scanTemporaryStorageBytes > 0)
                cudaErrCheck(cudaMalloc(&scanTemporaryStorage, scanTemporaryStorageBytes));

            cudaErrCheck(cudaMallocManaged(&obstacleParticlePairCount, sizeof(unsigned int)));
            cudaErrCheck(cudaStreamCreate(&obstacleStream));
            cudaErrCheck(cudaStreamCreate(&particleStream));
        }

        ~NeighborListLinkedCellState()
        {
            if(scanTemporaryStorage != nullptr)
                cudaErrCheck(cudaFree(scanTemporaryStorage));
            if(obstacleParticlePairCount != nullptr)
                cudaErrCheck(cudaFree(obstacleParticlePairCount));
            if(obstacleStream != nullptr)
                cudaErrCheck(cudaStreamDestroy(obstacleStream));
            if(particleStream != nullptr)
                cudaErrCheck(cudaStreamDestroy(particleStream));
        }
    };

    void NeighborListLinkedCellStateDeleter::operator()(NeighborListLinkedCellState* state) const
    {
        delete state;
    }

    NeighborListLinkedCellStatePtr createNeighborListLinkedCellState(unsigned int numParticles)
    {
        return NeighborListLinkedCellStatePtr(new NeighborListLinkedCellState(numParticles));
    }

    template <typename T>
    void updateLinkedCellNeighborList(NeighborListLinkedCellState&        state,
                                      NeighborListLinkedCellUpdateView<T> view)
    {
        using GP = GrainsParameters<T>;
        auto& CD = GP::m_collisionDetection;
        auto& LC = CD.linkedCellParameters;

        *view.pairCount = 0;

        const auto neighborCellsList = view.linkedCell->getCellNeighborsList();
        const auto cellParticleIDs   = view.linkedCell->getCellParticleIDs();
        const auto cellPrefixSums    = view.linkedCell->getCellPrefixSums();
        const auto numCells          = view.linkedCell->getNumCells();

        if(view.numObstacles > 0)
        {
            // Reset on the host before launch: block-local synchronization cannot order other
            // blocks before their atomic increments.
            *state.obstacleParticlePairCount = 0;
            generateObstacleParticlePairs_Device<<<view.numObstacles,
                                                   64,
                                                   0,
                                                   state.obstacleStream>>>(
                view.linkedCell->getObstacleIDs(),
                view.linkedCell->getObstacleCellIDs(),
                cellParticleIDs,
                cellPrefixSums,
                view.linkedCell->getMaxCellsPerObstacle(),
                view.numObstacles,
                view.numParticles,
                numCells,
                view.pairList->getData(),
                state.obstacleParticlePairCount);
            cudaStreamSynchronize(state.obstacleStream);
        }
        else
        {
            *state.obstacleParticlePairCount = 0;
        }

        unsigned int numBlocks;
        unsigned int numThreads;
        computeOptimalThreadsAndBlocks(view.numParticles,
                                       GrainsParameters<T>::m_GPU,
                                       numBlocks,
                                       numThreads);

        state.numNeighbors.fill(0u, state.particleStream);
        if(LC.type == LinkedCellType::ATOMICFIXED)
        {
            auto* linkedCellAtomicFixed = static_cast<LinkedCell_AtomicFixed<T>*>(view.linkedCell);
            countNeighbors_AtomicFixed_Device<<<numBlocks, numThreads, 0, state.particleStream>>>(
                neighborCellsList,
                linkedCellAtomicFixed->getCellParticleIDs(),
                linkedCellAtomicFixed->getAtomicPackBuffer(),
                linkedCellAtomicFixed->getNumParticlesPerCell(),
                linkedCellAtomicFixed->getMaxParticlesPerCell(),
                view.numParticles,
                view.numObstacles,
                numCells,
                state.numNeighbors.getData());
        }
        else
        {
            countNeighbors_Device<<<numBlocks, numThreads, 0, state.particleStream>>>(
                neighborCellsList,
                cellParticleIDs,
                cellPrefixSums,
                view.numParticles,
                view.numObstacles,
                numCells,
                state.numNeighbors.getData());
        }

        detail::exclusiveSumAsync(state.numNeighbors.getData(),
                                  state.numNeighborsPrefixSums.getData(),
                                  view.numParticles,
                                  state.scanTemporaryStorage,
                                  state.scanTemporaryStorageBytes,
                                  state.particleStream);

        unsigned int lastPrefixSum;
        unsigned int lastNeighborCount;
        cudaMemcpyAsync(&lastPrefixSum,
                        &state.numNeighborsPrefixSums.getData()[view.numParticles - 1],
                        sizeof(unsigned int),
                        cudaMemcpyDeviceToHost,
                        state.particleStream);
        cudaMemcpyAsync(&lastNeighborCount,
                        &state.numNeighbors.getData()[view.numParticles - 1],
                        sizeof(unsigned int),
                        cudaMemcpyDeviceToHost,
                        state.particleStream);

        cudaStreamSynchronize(state.obstacleStream);
        cudaStreamSynchronize(state.particleStream);

        unsigned int totalPairs = lastPrefixSum + lastNeighborCount;
        if(view.numObstacles > 0)
            totalPairs += *state.obstacleParticlePairCount;

        *view.pairCount = totalPairs;
        if(totalPairs != view.pairList->getSize())
            view.pairList->resize(totalPairs);

        if(LC.type == LinkedCellType::ATOMICFIXED)
        {
            auto* linkedCellAtomicFixed = static_cast<LinkedCell_AtomicFixed<T>*>(view.linkedCell);
            updateNeighborList_LC_AtomicFixed_Device<<<numBlocks,
                                                       numThreads,
                                                       0,
                                                       state.particleStream>>>(
                neighborCellsList,
                linkedCellAtomicFixed->getCellParticleIDs(),
                linkedCellAtomicFixed->getAtomicPackBuffer(),
                linkedCellAtomicFixed->getNumParticlesPerCell(),
                state.numNeighborsPrefixSums.getData(),
                linkedCellAtomicFixed->getMaxParticlesPerCell(),
                *state.obstacleParticlePairCount,
                view.numObstacles,
                view.numParticles,
                numCells,
                view.pairList->getData());
        }
        else
        {
            updateNeighborList_LC_Device<<<numBlocks, numThreads, 0, state.particleStream>>>(
                neighborCellsList,
                cellParticleIDs,
                cellPrefixSums,
                state.numNeighborsPrefixSums.getData(),
                *state.obstacleParticlePairCount,
                view.numObstacles,
                view.numParticles,
                numCells,
                view.pairList->getData());
        }

        cudaStreamSynchronize(state.particleStream);
    }

    template void updateLinkedCellNeighborList<float>(NeighborListLinkedCellState&,
                                                      NeighborListLinkedCellUpdateView<float>);
    template void updateLinkedCellNeighborList<double>(NeighborListLinkedCellState&,
                                                       NeighborListLinkedCellUpdateView<double>);
}
