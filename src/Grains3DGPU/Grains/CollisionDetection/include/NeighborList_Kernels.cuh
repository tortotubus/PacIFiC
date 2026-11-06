#ifndef NEIGHBOR_LIST_KERNELS_CUH
#define NEIGHBOR_LIST_KERNELS_CUH

#include "Basic.hh"

__GLOBAL__ void updateNeighborList_Nsq_Device(uint       nObstacles,
                                              uint       nParticles,
                                              IndexPair* pairList);

__GLOBAL__ void generateObstacleParticlePairs_Device(const IndexPair* obstacleIDs,
                                                      const uint*      obstacleCellIDs,
                                                      const uint64_t*  cellParticleIDs,
                                                      const uint*      cellPrefixSums,
                                                      uint             maxCellsPerObstacle,
                                                      uint             numObstacles,
                                                      uint             numParticles,
                                                      uint             numCells,
                                                      IndexPair*       pairList,
                                                      uint*            pairCount);

__GLOBAL__ void countNeighbors_Device(const uint*     cellNeighborsList,
                                      const uint64_t* cellParticleIDs,
                                      const uint*     cellPrefixSums,
                                      uint            numParticles,
                                      uint            numObstacles,
                                      uint            numCells,
                                      uint*           neighborCounts);

__GLOBAL__ void updateNeighborList_LC_Device(const uint*     cellNeighborsList,
                                             const uint64_t* cellParticleIDs,
                                             const uint*     cellPrefixSums,
                                             const uint*     numNeighborsPrefixSums,
                                             uint            pairListOffset,
                                             uint            numObstacles,
                                             uint            numParticles,
                                             uint            numCells,
                                             IndexPair*      pairList);

__GLOBAL__ void countNeighbors_AtomicFixed_Device(const uint*     cellNeighborsList,
                                                  const uint64_t* cellParticleIDsFixed,
                                                  const uint64_t* particleCellIDsSeq,
                                                  const uint*     numParticlesPerCell,
                                                  uint            maxParticlesPerCell,
                                                  uint            numParticles,
                                                  uint            numObstacles,
                                                  uint            numCells,
                                                  uint*           neighborCounts);

__GLOBAL__ void updateNeighborList_LC_AtomicFixed_Device(
    const uint*     cellNeighborsList,
    const uint64_t* cellParticleIDsFixed,
    const uint64_t* particleCellIDsSeq,
    const uint*     numParticlesPerCell,
    const uint*     numNeighborsPrefixSums,
    uint            maxParticlesPerCell,
    uint            pairListOffset,
    uint            numObstacles,
    uint            numParticles,
    uint            numCells,
    IndexPair*      pairList);

#endif
