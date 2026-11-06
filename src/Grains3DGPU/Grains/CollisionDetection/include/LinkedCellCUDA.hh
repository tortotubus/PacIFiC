#ifndef GRAINS_LINKED_CELL_CUDA_HH
#define GRAINS_LINKED_CELL_CUDA_HH

#include <cstddef>
#include <cstdint>
#include <memory>

#include "CellsFwd.hh"
#include "GrainsMemBuffer.hh"
#include "GrainsTypes.hh"

template <typename T>
class Quaternion;
template <typename T>
class RigidBody;
template <typename T>
class Vector3;

namespace grains::cuda
{
    namespace detail
    {
        std::size_t exclusiveSumStorageBytes(const unsigned int* input,
                                             unsigned int*       output,
                                             unsigned int        count);
        void        exclusiveSum(const unsigned int* input,
                                 unsigned int*       output,
                                 unsigned int        count,
                                 void*               scratch,
                                 std::size_t         scratchBytes);
        void        exclusiveSumAsync(const unsigned int* input,
                                      unsigned int*       output,
                                      unsigned int        count,
                                      void*               scratch,
                                      std::size_t         scratchBytes,
                                      void*               stream);

    }

    /** Host-side inputs used by the CUDA SortBased linked-cell update. */
    template <typename T>
    struct SortBasedUpdateView
    {
        Cells<T>**                                      cells         = nullptr;
        const Vector3<T>*                               positions     = nullptr;
        Vector3<T>*                                     oldPositions  = nullptr;
        GrainsMemBuffer<unsigned int, MemType::DEVICE>* neighborCells = nullptr;
        unsigned int                                    numCells      = 0;
        unsigned int                                    numObstacles  = 0;
        unsigned int                                    numParticles  = 0;
        bool                                            adaptive      = false;
        T                                               cellSize      = T(0);
    };

    /** CUDA-owned streams, event, workspace, and output buffers for SortBased linked cells. */
    struct SortBasedState;
    struct SortBasedStateDeleter
    {
        void operator()(SortBasedState* state) const;
    };
    using SortBasedStatePtr = std::unique_ptr<SortBasedState, SortBasedStateDeleter>;

    SortBasedStatePtr    createSortBasedState(unsigned int numParticles, unsigned int numCells);
    const std::uint64_t* sortBasedCellParticleIDs(const SortBasedState& state);
    const unsigned int*  sortBasedCellPrefixSums(const SortBasedState& state);

    /** Prepare hashes and optional adaptive grid resizing; returns the current cell count. */
    template <typename T>
    unsigned int prepareSortBasedUpdate(SortBasedState& state, SortBasedUpdateView<T> view);

    /** Sort packed particle IDs and rebuild the cell-prefix table. */
    template <typename T>
    void finalizeSortBasedUpdate(SortBasedState& state, unsigned int numParticles);

    /** CUDA-owned resources used to reduce particle displacements. */
    template <typename T>
    struct LinkedCellReductionState;

    struct LinkedCellReductionStateDeleter
    {
        void operator()(LinkedCellReductionState<float>* state) const;
        void operator()(LinkedCellReductionState<double>* state) const;
    };

    template <typename T>
    using LinkedCellReductionStatePtr
        = std::unique_ptr<LinkedCellReductionState<T>, LinkedCellReductionStateDeleter>;

    /** Allocate the CUDA stream and reduction workspace for at most capacity positions. */
    template <typename T>
    LinkedCellReductionStatePtr<T> createLinkedCellReductionState(std::size_t       capacity,
                                                                  const Vector3<T>* oldPositions,
                                                                  const Vector3<T>* positions,
                                                                  T*                output);

    /** Compute the maximum squared displacement and wait for the result. */
    template <typename T>
    void computeMaxDisplacementSquared(LinkedCellReductionState<T>& state,
                                       const Vector3<T>*            oldPositions,
                                       const Vector3<T>*            positions,
                                       std::size_t                  count,
                                       T*                           output);

    /** Launch the maximum squared displacement reduction without waiting for it. */
    template <typename T>
    void launchMaxDisplacementSquared(LinkedCellReductionState<T>& state,
                                      const Vector3<T>*            oldPositions,
                                      const Vector3<T>*            positions,
                                      std::size_t                  count,
                                      T*                           output);

    /** Wait for work submitted to the linked-cell reduction stream. */
    template <typename T>
    void synchronizeLinkedCellReduction(LinkedCellReductionState<T>& state);

    /** Launch the CUDA neighbor-cell generation kernel and wait for completion. */
    template <typename T>
    void generateNeighborCells(const Cells<T>* const* cells,
                               unsigned int           numCells,
                               unsigned int*          neighborCells);

    /** Launch the CUDA obstacle-linking kernel and wait for completion. */
    template <typename T>
    void linkObstacles(const RigidBody<T>* const* rigidBodies,
                       const Vector3<T>*          positions,
                       const Quaternion<T>*       quaternions,
                       const Cells<T>* const*     cells,
                       unsigned int               numObstacles,
                       unsigned int               maxCellsPerObstacle,
                       IndexPair*                 obstacleIDs,
                       unsigned int*              obstacleCellIDs);
}

#endif
