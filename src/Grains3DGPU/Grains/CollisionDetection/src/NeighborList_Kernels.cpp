#include "NeighborList_Kernels.hh"

// -------------------------------------------------------------------------------------------------
// Updates the neighbor list on host using an O(n^2) algorithm
__HOST__ void
    updateNeighborList_Nsq_Host(const uint nObstacles, const uint nParticles, IndexPair* pairList)
{
    for(uint i = 0; i < nObstacles; ++i)
        for(uint j = 0; j < nParticles; ++j)
            pairList[nParticles * i + j] = makeIndexPair(i, nObstacles + j);

    // Offset for p-p interactions.
    uint offset = nObstacles * nParticles;
    for(uint i = 0; i < nParticles; ++i)
        for(uint j = i + 1; j < nParticles; ++j)
            pairList[offset + i + j * (j - 1) / 2] = makeIndexPair(nObstacles + i, nObstacles + j);
}

// -------------------------------------------------------------------------------------------------
// Updates the neighbor list on host using a linked cell approach
__HOST__ void updateNeighborList_LC_Host(const uint*                         cellNeighborsList,
                                         const IndexPair*                    obstacleIDs,
                                         const uint*                         obstacleCellIDs,
                                         const uint*                         particleIDs,
                                         const uint*                         cellIDs,
                                         const std::vector<std::list<uint>>& cellParticles,
                                         const uint                          maxCellsPerObstacle,
                                         const uint                          numObstacles,
                                         const uint                          numParticles,
                                         GrainsMemBuffer<IndexPair, MemType::HOST>& pairList)
{
    constexpr uint NUM_NEIGHBOR_CELLS = 27;  // Number of neighboring cells

    // Set pairList size to zero
    pairList.setSize(0);

    // FIRST PASS: Loop over all obstacles
    for(uint i = 0; i < numObstacles; ++i)
    {
        const uint offset = i * maxCellsPerObstacle;

        const uint obstacleIndex      = obstacleIDs[i].x;
        const uint numCellsToTraverse = obstacleIDs[i].y;

        for(uint c = 0; c < numCellsToTraverse; ++c)
        {
            const uint  cell                = obstacleCellIDs[offset + c];
            const auto& targetCellParticles = cellParticles[cell];

            // Check against all particles in the target cell
            for(uint particleID : targetCellParticles)
                pairList.push_back(makeIndexPair(obstacleIndex, particleID));
        }
    }

    // SECOND PASS: Cell-centric approach for particle-particle pairs
    for(uint c = 0; c < cellParticles.size(); ++c)
    {
        if(cellParticles[c].empty())
            continue;  // Skip empty cells

        // Get neighbor cells for this cell (includes own cell)
        const uint* neighborCells = &cellNeighborsList[NUM_NEIGHBOR_CELLS * c];

        // Check interactions with neighboring cells
        for(uint cc = 0; cc < NUM_NEIGHBOR_CELLS; ++cc)
        {
            uint targetCell = neighborCells[cc];
            if(targetCell == UINT_MAX || targetCell < c)
                continue;  // Skip invalid cells and cells with lower indices

            const auto& neighborCellParticles = cellParticles[targetCell];
            if(neighborCellParticles.empty())
                continue;  // Skip empty target cells

            // Process particle pairs between cells
            for(uint primaryParticle : cellParticles[c])
            {
                for(uint otherParticle : neighborCellParticles)
                {
                    if(targetCell == c)
                    {
                        // Same-cell: avoid duplicates with ordering check
                        if(primaryParticle < otherParticle)
                            pairList.push_back(makeIndexPair(primaryParticle, otherParticle));
                    }
                    else
                    {
                        // Cross-cell: always emit -- uniqueness is guaranteed by targetCell > c
                        const uint lo
                            = primaryParticle < otherParticle ? primaryParticle : otherParticle;
                        const uint hi
                            = primaryParticle < otherParticle ? otherParticle : primaryParticle;
                        pairList.push_back(makeIndexPair(lo, hi));
                    }
                }
            }
        }
    }
}
