#ifndef GRAINS_NEIGHBOR_LIST_LINKED_CELL_CUDA_HH
#define GRAINS_NEIGHBOR_LIST_LINKED_CELL_CUDA_HH

#include <memory>

#include "GrainsMemBuffer.hh"
#include "GrainsTypes.hh"

template <typename T, MemType M>
class LinkedCell;

namespace grains::cuda
{
    /** CUDA-owned buffers, scan workspace, counter, and streams for linked-cell neighbor lists. */
    struct NeighborListLinkedCellState;
    struct NeighborListLinkedCellStateDeleter
    {
        void operator()(NeighborListLinkedCellState* state) const;
    };
    using NeighborListLinkedCellStatePtr
        = std::unique_ptr<NeighborListLinkedCellState, NeighborListLinkedCellStateDeleter>;

    NeighborListLinkedCellStatePtr createNeighborListLinkedCellState(unsigned int numParticles);

    /** Shared inputs and outputs needed by the CUDA neighbor-list update. */
    template <typename T>
    struct NeighborListLinkedCellUpdateView
    {
        LinkedCell<T, MemType::DEVICE>*              linkedCell   = nullptr;
        GrainsMemBuffer<IndexPair, MemType::DEVICE>* pairList     = nullptr;
        unsigned int*                                pairCount    = nullptr;
        unsigned int                                 numObstacles = 0;
        unsigned int                                 numParticles = 0;
    };

    /** Rebuild the device pair list from the current linked-cell data. */
    template <typename T>
    void updateLinkedCellNeighborList(NeighborListLinkedCellState&        state,
                                      NeighborListLinkedCellUpdateView<T> view);
}

#endif
