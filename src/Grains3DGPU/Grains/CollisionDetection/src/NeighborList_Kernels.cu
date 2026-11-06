#include "GrainsParameters.hh"
#include "GrainsPlatform.hh"
#include "GrainsUtils.hh"
#include "NeighborList_Kernels.cuh"
#include "NeighborList_Kernels.hh"

// -------------------------------------------------------------------------------------------------
// Updates the neighbor list on device using an O(n^2) algorithm

__GLOBAL__ void
    updateNeighborList_Nsq_Device(const uint nObstacles, const uint nParticles, IndexPair* pairList)
{
    uint tID = blockIdx.x * blockDim.x + threadIdx.x;
    if(tID < nObstacles)
    {
        // Obstacle to particle pairs
        for(uint j = 0; j < nParticles; ++j)
            pairList[nParticles * tID + j] = makeIndexPair(tID, nObstacles + j);
    }
    else if(tID < nObstacles + nParticles)
    {
        // offset
        const uint offset = nObstacles * nParticles;
        // adjust tID to start from 0 for particles
        tID -= nObstacles;
        // Particle to obstacle pairs
        for(uint j = tID + 1; j < nParticles; ++j)
            pairList[offset + tID + j * (j - 1) / 2]
                = makeIndexPair(nObstacles + tID, nObstacles + j);
    }
    else
        return;
}

namespace grains::cuda
{
    template <typename T>
    void updateNeighborListNsq(uint nObstacles, uint nParticles, IndexPair* pairList)
    {
        uint numBlocks;
        uint numThreads;
        computeOptimalThreadsAndBlocks(nObstacles + nParticles,
                                       GrainsParameters<T>::m_GPU,
                                       numBlocks,
                                       numThreads);
        updateNeighborList_Nsq_Device<<<numBlocks, numThreads>>>(nObstacles, nParticles, pairList);
        cudaDeviceSynchronize();
    }

    template void updateNeighborListNsq<float>(uint, uint, IndexPair*);
    template void updateNeighborListNsq<double>(uint, uint, IndexPair*);
}

// -------------------------------------------------------------------------------------------------
// Generate obstacle-particle pairs on device

__GLOBAL__ void generateObstacleParticlePairs_Device(const IndexPair* obstacleIDs,
                                                     const uint*      obstacleCellIDs,
                                                     const uint64_t*  packedCellParticleIDs,
                                                     const uint*      cellPrefixSums,
                                                     const uint       maxCellsPerObstacle,
                                                     const uint       numObstacles,
                                                     const uint       numParticles,
                                                     const uint       numCells,
                                                     IndexPair*       pairList,
                                                     uint*            pairCount)
{
    uint obstacleIdx = blockIdx.x;
    if(obstacleIdx >= numObstacles)
        return;

    const uint offset             = obstacleIdx * maxCellsPerObstacle;
    const uint obstacleIndex      = obstacleIDs[obstacleIdx].x;
    const uint numCellsToTraverse = obstacleIDs[obstacleIdx].y;

    // Each thread handles one cell for this obstacle
    for(uint c = threadIdx.x; c < numCellsToTraverse; c += blockDim.x)
    {
        const uint cell      = obstacleCellIDs[offset + c];
        const uint cellStart = cellPrefixSums[cell];

        if(cellStart == UINT_MAX)
            continue;  // Empty cell

        // Find cell end (handles both variants: with/without UINT_MAX sentinels)
        uint cellEnd;
        uint k = cell;
        do
        {
            ++k;
            cellEnd = (k < numCells) ? cellPrefixSums[k] : numParticles;
        } while(cellEnd == UINT_MAX && k < numCells);

        // Add pairs for all particles in this cell
        for(uint p = cellStart; p < cellEnd; ++p)
        {
            // Unpack particleID from lower 32 bits
            uint64_t packed = packedCellParticleIDs[p];
            // Skip empty slots: AtomicFixed layout fills unused entries with UINT64_MAX.
            // Without this check, all maxParticlesPerCell slots are written as pairs
            // (most with particleID = UINT_MAX), causing a massive OOB overwrite of the
            // pair list buffer.
            if(packed == UINT64_MAX)
                continue;
            uint particleID       = static_cast<uint>(packed & 0xFFFFFFFF);
            uint globalIndex      = atomicAdd(pairCount, 1);
            pairList[globalIndex] = makeIndexPair(obstacleIndex, particleID);
        }
    }
}

// -------------------------------------------------------------------------------------------------
// Count neighbors per particle using linked cells with packed keys

__GLOBAL__ void countNeighbors_Device(const uint*     cellNeighborsList,
                                      const uint64_t* cellParticleIDs,
                                      const uint*     cellPrefixSums,
                                      const uint      numParticles,
                                      const uint      numObstacles,
                                      const uint      numCells,
                                      uint*           neighborCounts)
{
    constexpr uint NUM_NEIGHBOR_CELLS = 27;

    uint tID = blockIdx.x * blockDim.x + threadIdx.x;
    if(tID >= numParticles)
        return;

    // Unpack cellID and particleID from uint64
    const uint64_t packed = cellParticleIDs[tID];
    const uint     cell   = (uint)(packed >> 32);
    const uint     i      = (uint)(packed & 0xFFFFFFFF);  // global particle ID

    if(cell == UINT_MAX)
        return;  // neighborCounts pre-filled with 0; no write needed

    const uint* neighborCells  = &cellNeighborsList[NUM_NEIGHBOR_CELLS * cell];
    uint        totalNeighbors = 0;

    // Loop over all neighboring cells
    for(uint cID = 0; cID < NUM_NEIGHBOR_CELLS; ++cID)
    {
        uint c = neighborCells[cID];

        if(c == UINT_MAX || c < cell)
            continue;

        // Calculate number of particles in cell from cellPrefixSums
        uint cellStart = cellPrefixSums[c];
        if(cellStart == UINT_MAX)
            continue;

        // Find end of cell
        uint k = c + 1;
        while(k < numCells && cellPrefixSums[k] == UINT_MAX)
            ++k;
        const uint cellEnd = (k < numCells) ? cellPrefixSums[k] : numParticles;

        // Count particles in this neighbor cell.
        // For same-cell pairs (c == cell): use i < j to avoid duplicates.
        // For cross-cell pairs (c > cell): the c >= cell condition already ensures
        // only this thread processes the (cell, c) cell pair, so emit all pairs.
        for(uint p = cellStart; p < cellEnd; ++p)
        {
            const uint64_t neighbor = cellParticleIDs[p];
            if(neighbor == UINT64_MAX)
                break;
            const uint j = (uint)(neighbor & 0xFFFFFFFF);

            if(c == cell && j <= i)
                continue;

            totalNeighbors++;
        }
    }

    neighborCounts[i - numObstacles] = totalNeighbors;
}

// -------------------------------------------------------------------------------------------------
// Updates the neighbor list on device

__GLOBAL__ void updateNeighborList_LC_Device(const uint*     cellNeighborsList,
                                             const uint64_t* cellParticleIDs,
                                             const uint*     cellPrefixSums,
                                             const uint*     numNeighborsPrefixSums,
                                             const uint      pairListOffset,
                                             const uint      numObstacles,
                                             const uint      numParticles,
                                             const uint      numCells,
                                             IndexPair*      pairList)
{
    constexpr uint NUM_NEIGHBOR_CELLS = 27;

    uint tID = blockIdx.x * blockDim.x + threadIdx.x;
    if(tID >= numParticles)
        return;

    // Unpack cellID and particleID from uint64
    const uint64_t packed = cellParticleIDs[tID];
    const uint     cell   = (uint)(packed >> 32);
    const uint     i      = (uint)(packed & 0xFFFFFFFF);

    if(cell == UINT_MAX || cell >= numCells)
        return;

    const uint* neighborCells = &cellNeighborsList[NUM_NEIGHBOR_CELLS * cell];
    uint        insertIndex   = pairListOffset + numNeighborsPrefixSums[i - numObstacles];

    // Loop over all neighboring cells
    for(uint cID = 0; cID < NUM_NEIGHBOR_CELLS; ++cID)
    {
        uint c = neighborCells[cID];

        if(c == UINT_MAX || c < cell)
            continue;

        uint cellStart = cellPrefixSums[c];
        if(cellStart == UINT_MAX)
            continue;

        // Find end of cell
        uint k = c + 1;
        while(k < numCells && cellPrefixSums[k] == UINT_MAX)
            ++k;
        const uint cellEnd = (k < numCells) ? cellPrefixSums[k] : numParticles;

        // Loop through particles in the neighbor cell.
        // For same-cell pairs (c == cell): use i < j to avoid duplicates.
        // For cross-cell pairs (c > cell): the c >= cell condition already ensures
        // only this thread processes the (cell, c) cell pair, so emit all pairs.
        for(uint p = cellStart; p < cellEnd; ++p)
        {
            // Unpack particle ID from cellParticleIDs
            const uint64_t neighbor = cellParticleIDs[p];
            if(neighbor == UINT64_MAX)
                break;
            const uint j = (uint)(neighbor & 0xFFFFFFFF);

            if(c == cell && i >= j)
                continue;

            // Always write (min,max) to keep pairs canonically ordered.
            pairList[insertIndex++] = (i < j) ? makeIndexPair(i, j) : makeIndexPair(j, i);
        }
    }
}

// -------------------------------------------------------------------------------------------------
// Count neighbors per particle for AtomicFixed using direct count lookup

__GLOBAL__ void countNeighbors_AtomicFixed_Device(const uint*     cellNeighborsList,
                                                  const uint64_t* cellParticleIDsFixed,
                                                  const uint64_t* particleCellIDsSeq,
                                                  const uint*     numParticlesPerCell,
                                                  const uint      maxParticlesPerCell,
                                                  const uint      numParticles,
                                                  const uint      numObstacles,
                                                  const uint      numCells,
                                                  uint*           neighborCounts)
{
    constexpr uint NUM_NEIGHBOR_CELLS = 27;

    uint tID = blockIdx.x * blockDim.x + threadIdx.x;
    if(tID >= numParticles)
        return;

    // Use sequential buffer: particleCellIDsSeq[tID] = (cellID<<32 | particleID)
    // cellParticleIDsFixed uses 2D layout [cellID*maxPerCell+slot] for neighbour traversal.
    const uint64_t packed = particleCellIDsSeq[tID];
    const uint     cell   = (uint)(packed >> 32);
    const uint     i      = (uint)(packed & 0xFFFFFFFF);

    if(cell == UINT_MAX || cell >= numCells)
        return;

    const uint* neighborCells = &cellNeighborsList[NUM_NEIGHBOR_CELLS * cell];
    uint        count         = 0;

    // Loop over all neighboring cells
    for(uint cID = 0; cID < NUM_NEIGHBOR_CELLS; ++cID)
    {
        uint c = neighborCells[cID];

        if(c == UINT_MAX || c < cell)
            continue;

        // Compute cellStart for fixed layout: cellID * maxParticlesPerCell
        uint cellStart = c * maxParticlesPerCell;
        uint cellEnd   = cellStart + numParticlesPerCell[c];

        // Count particles in the neighbor cell.
        // For same-cell pairs (c == cell): use i < j to avoid duplicates.
        // For cross-cell pairs (c > cell): the c >= cell condition already ensures
        // only this thread processes the (cell, c) cell pair, so emit all pairs.
        for(uint p = cellStart; p < cellEnd; ++p)
        {
            const uint j = (uint)(cellParticleIDsFixed[p] & 0xFFFFFFFF);

            if(c == cell && i >= j)
                continue;

            count++;
        }
    }

    neighborCounts[i - numObstacles] = count;
}

// -------------------------------------------------------------------------------------------------
// Updates the neighbor list for AtomicFixed using direct count lookup

__GLOBAL__ void updateNeighborList_LC_AtomicFixed_Device(const uint*     cellNeighborsList,
                                                         const uint64_t* cellParticleIDsFixed,
                                                         const uint64_t* particleCellIDsSeq,
                                                         const uint*     numParticlesPerCell,
                                                         const uint*     numNeighborsPrefixSums,
                                                         const uint      maxParticlesPerCell,
                                                         const uint      pairListOffset,
                                                         const uint      numObstacles,
                                                         const uint      numParticles,
                                                         const uint      numCells,
                                                         IndexPair*      pairList)
{
    constexpr uint NUM_NEIGHBOR_CELLS = 27;

    uint tID = blockIdx.x * blockDim.x + threadIdx.x;
    if(tID >= numParticles)
        return;

    // Use sequential buffer: particleCellIDsSeq[tID] = (cellID<<32 | particleID)
    // cellParticleIDsFixed uses 2D layout [cellID*maxPerCell+slot] for neighbour traversal.
    const uint64_t packed = particleCellIDsSeq[tID];
    const uint     cell   = (uint)(packed >> 32);
    const uint     i      = (uint)(packed & 0xFFFFFFFF);

    if(cell == UINT_MAX || cell >= numCells)
        return;

    const uint* neighborCells = &cellNeighborsList[NUM_NEIGHBOR_CELLS * cell];
    uint        insertIndex   = pairListOffset + numNeighborsPrefixSums[i - numObstacles];

    // Loop over all neighboring cells
    for(uint cID = 0; cID < NUM_NEIGHBOR_CELLS; ++cID)
    {
        uint c = neighborCells[cID];

        if(c == UINT_MAX || c < cell)
            continue;

        // Compute cellStart for fixed layout: cellID * maxParticlesPerCell
        uint cellStart = c * maxParticlesPerCell;
        uint cellEnd   = cellStart + numParticlesPerCell[c];

        // Loop through particles in the neighbor cell.
        // For same-cell pairs (c == cell): use i < j to avoid duplicates.
        // For cross-cell pairs (c > cell): the c >= cell condition already ensures
        // only this thread processes the (cell, c) cell pair, so emit all pairs.
        for(uint p = cellStart; p < cellEnd; ++p)
        {
            const uint j = (uint)(cellParticleIDsFixed[p] & 0xFFFFFFFF);

            if(c == cell && i >= j)
                continue;

            // Always write (min,max) to keep pairs canonically ordered.
            pairList[insertIndex++] = (i < j) ? makeIndexPair(i, j) : makeIndexPair(j, i);
        }
    }
}
