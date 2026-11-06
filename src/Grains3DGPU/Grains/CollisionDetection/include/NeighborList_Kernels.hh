#ifndef NEIGHBOR_LIST_KERNELS_HH
#define NEIGHBOR_LIST_KERNELS_HH

#include "Basic.hh"
#include "GrainsMemBuffer.hh"

// Host neighbor-list implementations and C++-parseable CUDA backend entry points. CUDA kernel
// declarations live in NeighborList_Kernels.cuh and are not part of the C++ interface.
__HOST__ void updateNeighborList_Nsq_Host(uint nObstacles, uint nParticles, IndexPair* pairList);

namespace grains::cuda
{
    template <typename T>
    void updateNeighborListNsq(uint nObstacles, uint nParticles, IndexPair* pairList);
}

__HOST__ void updateNeighborList_LC_Host(const uint*                         cellNeighborsList,
                                         const IndexPair*                    obstacleIDs,
                                         const uint*                         obstacleCellIDs,
                                         const uint*                         particleIDs,
                                         const uint*                         cellIDs,
                                         const std::vector<std::list<uint>>& cellParticles,
                                         uint                                maxCellsPerObstacle,
                                         uint                                numObstacles,
                                         uint                                numParticles,
                                         GrainsMemBuffer<IndexPair, MemType::HOST>& pairList);

#endif
