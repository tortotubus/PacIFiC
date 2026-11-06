#ifndef _LINKEDCELL_HH_
#define _LINKEDCELL_HH_

#include "GrainsPlatform.hh"

#include "Cells.hh"
#include "CellsFactory.hh"
#include "GrainsMemBuffer.hh"
#include "GrainsUtils.hh"
#include "LinkedCellCUDA.hh"
#include "QuaternionMath.hh"
#include "VectorMath.hh"

namespace grains::detail
{
    template <typename T, MemType M>
    struct LinkedCellReductionStorage
    {
    };

    template <typename T>
    struct LinkedCellReductionStorage<T, MemType::DEVICE>
    {
        cuda::LinkedCellReductionStatePtr<T> state;
    };
}

// =================================================================================================
/** @brief The class LinkedCell.

    This class provides functionalities to to manage linked cells for
    collision detection in the simulation. It is essentially a wrapper around
    the LinkedCell class, providing methods to create and update the neighbor
    list based on the linked cells. This wrapper is designed to work only on
    host.
    Note that number of obstacles is fixed after construction. This is a hard
    constraint since we allocate buffers based on this number.

    @author A.Yazdani - 2025 - Construction */
// =================================================================================================
template <typename T, MemType M>
class LinkedCell
{
    static_assert(M == MemType::HOST || M == MemType::DEVICE,
                  "LinkedCell only supports MemType::HOST or MemType::DEVICE");

protected:
    /** @name Parameters */
    //@{
    /** \brief Non-owning pointer to the rigid bodies buffer (stable address). We assume that this
        buffer remains valid during the lifetime of this object. */
    const GrainsMemBuffer<RigidBody<T>*, M>* m_rb = nullptr;
    /** \brief Non-owning pointer to positions buffer */
    const GrainsMemBuffer<Vector3<T>, M>* m_positions = nullptr;
    /** \brief Non-owning pointer to quaternions buffer */
    const GrainsMemBuffer<Quaternion<T>, M>* m_quaternions = nullptr;
    /** \brief Particles position in the last update */
    GrainsMemBuffer<Vector3<T>, M> m_oldPosition;
    /** \brief Cells object. We allocate a buffer for later if we want to work with multiple cells
        objects. */
    GrainsMemBuffer<Cells<T>*, M> m_cells;
    /** \brief Buffer to store neighbor cell IDs */
    GrainsMemBuffer<uint, M> m_neighborCells;
    /** \brief Buffer of obstacle IDs and the number of cells that have to be checked for a possible
        contact with a particle. This is essentially the number of cells each obstacle occupies +
        one-ring. */
    GrainsMemBuffer<IndexPair, M> m_obstacleID;
    /** \brief Buffer of the cell IDs that that have to be checked for a possible contact with an
        obstacle. */
    GrainsMemBuffer<uint, M> m_obstacleCellID;
    /** \brief CUDA-owned state for maximum-displacement reduction. */
    grains::detail::LinkedCellReductionStorage<T, M> m_reductionBackend;
    /** \brief Maximum displacement of particles since the last update. Note that we store the
        squared value.  */
    GrainsMemBuffer<T, M == MemType::HOST ? MemType::HOST : MemType::MANAGED>
        m_maxDisplacementSquared;
    /** \brief Cell size. This is the minimum possible size for the cells. Skin thickness will be
        added to this value. */
    T m_cellSizeWithoutSkin = 0;
    /** \brief Skin thickness */
    T m_skinThickness = 0;
    /** \brief Update frequency */
    uint m_updateFrequency = 0;
    /** \brief Number of iterations since the last update */
    uint m_numIterationsSinceLastUpdate = 0;
    /** \brief Number of obstacles */
    uint m_numObstacles = 0;
    /** \brief Maximum number of cells an obstacle can occupy + 1-ring */
    uint m_maxCellsPerObstacle = 0;
    /** \brief Number of particles */
    uint m_numParticles = 0;
    /** \brief Number of cells in the grid */
    uint m_numCells = 0;
    /** \brief Flag to indicate if adaptive skin is used */
    bool m_useAdaptiveSkin = false;
    /** \brief Whether a valid async displacement result is available */
    bool m_hasAsyncResult = false;  //@}

public:
    /** @name Constructors */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Default constructor */
    LinkedCell() = default;

    // ---------------------------------------------------------------------------------------------
    /** @brief Constructor with parameters
        @param rb Rigid body buffer
        @param positions Positions buffer
        @param quaternions Quaternions buffer
        @param linkedCellParameters Linked cell parameters
        @param nObstacles number of obstacles
        @param nParticles number of particles
        @param nCellsForEachObstacle number of cells for each obstacle */
    LinkedCell(const GrainsMemBuffer<RigidBody<T>*, M>* rb,
               const GrainsMemBuffer<Vector3<T>, M>&    positions,
               const GrainsMemBuffer<Quaternion<T>, M>& quaternions,
               const LinkedCellParameters<T>&           linkedCellParameters,
               const uint                               nObstacles,
               const uint                               nParticles)
        : m_oldPosition(nObstacles + nParticles)
        , m_obstacleID(nObstacles)
        , m_maxDisplacementSquared(1)
        , m_numObstacles(nObstacles)
        , m_numParticles(nParticles)
    {
        // Store non-owning pointer to rigid body buffer (must remain valid)
        m_rb          = rb;
        m_positions   = &positions;
        m_quaternions = &quaternions;

        // Sanity check
        GAssert(positions.getSize() == nObstacles + nParticles
                    && quaternions.getSize() == nObstacles + nParticles,
                "LinkedCell: positions or quaternions size does not match the number of bodies!");

        // Extract parameters
        const Vector3<T>& minCorner              = linkedCellParameters.minCorner;
        const Vector3<T>& maxCorner              = linkedCellParameters.maxCorner;
        const T           minCellSize            = linkedCellParameters.minCellSize;
        const T           cellSizeFactor         = linkedCellParameters.cellSizeFactor;
        const uint        maxNumCellsPerObstacle = linkedCellParameters.maxNumCellsPerObstacle;
        m_updateFrequency                        = linkedCellParameters.updateFrequency;
        m_useAdaptiveSkin                        = (m_updateFrequency > 0);

        // The minimum cell size should be at least twice the maximum radius
        // of particles. We multiply this value by a factor (>=1) to get the
        // final cell size.
        m_cellSizeWithoutSkin = minCellSize * cellSizeFactor;

        // Initialize the LinkedCell buffer
        // Note that we start with smallest size possible (largest number of cells) so the buffers
        // are allocated with the largest size possible.
        m_cells.initialize(1);
        m_numCells = CellsFactory<T>::template create<M>(minCorner,
                                                         maxCorner,
                                                         m_cellSizeWithoutSkin,
                                                         m_cells);

        // Initialize neighbor cells buffer
        m_neighborCells.initialize(m_numCells * 27);  // 26 neighbors + self
        m_neighborCells.fill(UINT_MAX);
        generateNeighborCells();

        // Compute initial skin thickness
        m_skinThickness = 0.1 * m_cellSizeWithoutSkin;

        // Adjust the maximum number of cells per obstacle
        m_maxCellsPerObstacle = std::min(maxNumCellsPerObstacle, m_numCells);
        m_obstacleCellID.initialize(m_maxCellsPerObstacle * nObstacles);

        // Initialize old positions with current positions for first displacement computation
        m_oldPosition.copyFrom(positions);

        // Setup obstacle
        linkObstacles();

        // Initialize the CUDA-owned reduction state for maximum displacement computation.
        if constexpr(M == MemType::DEVICE)
        {
            m_reductionBackend.state = grains::cuda::createLinkedCellReductionState<T>(
                nObstacles + nParticles,
                m_oldPosition.getData(),
                m_positions->getData(),
                m_maxDisplacementSquared.getData());
        }
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Destructor */
    virtual ~LinkedCell()
    {
        // Clean up the Cells objects using GrainsMemBuffer helper method
        // m_cells.freePointedObjects();

        // Set non-owning pointers to nullptr
        m_rb          = nullptr;
        m_positions   = nullptr;
        m_quaternions = nullptr;
    }
    //@}

    /** @name Get methods */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Gets linked cell list */
    Cells<T>* const* getLinkedCell() const
    {
        return m_cells.getData();
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Gets neighbor cells */
    const uint* getCellNeighborsList() const
    {
        return m_neighborCells.getData();
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Gets obstacle IDs */
    const IndexPair* getObstacleIDs() const
    {
        return m_obstacleID.getData();
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Gets obstacle cell IDs */
    const uint* getObstacleCellIDs() const
    {
        return m_obstacleCellID.getData();
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Gets cell size without skin thickness. */
    T getCellSizeWithoutSkin() const
    {
        return m_cellSizeWithoutSkin;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Gets skin thickness. */
    T getSkinThickness() const
    {
        return m_skinThickness;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Gets maximum displacement. */
    T getMaxDisplacement() const
    {
        return sqrt(m_maxDisplacementSquared);
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Gets number of iterations since last update. */
    uint getNumIterationsSinceLastUpdate() const
    {
        return m_numIterationsSinceLastUpdate;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Gets maximum number of cells an obstacle can occupy. */
    uint getMaxCellsPerObstacle() const
    {
        return m_maxCellsPerObstacle;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Gets number of cells. */
    uint getNumCells() const
    {
        return m_numCells;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Gets particle IDs (implementation-specific) */
    virtual const uint* getParticleIDs() const
    {
        GAbort("LinkedCell::getParticleIDs is not supported in this variant");
        return nullptr;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Gets cell IDs (implementation-specific) */
    virtual const uint* getCellIDs() const
    {
        GAbort("LinkedCell::getCellIDs is not supported in this variant");
        return nullptr;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Gets number of particles per cell (implementation-specific) */
    virtual const uint* getNumParticlesPerCell() const
    {
        GAbort("LinkedCell::getNumParticlesPerCell is not supported in this variant");
        return nullptr;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Gets packed cell-particle IDs (implementation-specific) */
    virtual const uint64_t* getCellParticleIDs() const
    {
        GAbort("LinkedCell::getCellParticleIDs is not supported in this variant");
        return nullptr;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Gets cell prefix sums (start indices for each cell in particle arrays) */
    virtual const uint* getCellPrefixSums() const
    {
        GAbort("LinkedCell::getCellPrefixSums is not supported in this variant");
        return nullptr;
    }
    //@}

    /** @name Set methods */
    //@{
    //@}

    /** @name Methods */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Checks if linked cell update is needed based on adaptive skin and displacement. */
    bool needsUpdate()
    {
        auto& SS = GrainsParameters<T>::m_simulationState;

        if(!m_useAdaptiveSkin)
            return true;  // Always update in non-adaptive mode

        // Increment before any early return so computeSkinThickness() never divides by zero
        ++m_numIterationsSinceLastUpdate;

        if(SS.particlesSorted)
            return true;  // Always update if particles were sorted

        if constexpr(M == MemType::HOST)
        {
            T maxDispSq = computeMaxDisplacement();
            return (T(4) * maxDispSq > m_skinThickness * m_skinThickness);
        }
        else
        {
            // 1-lag async: reduce was launched at end of prior updateLinkedCells(); sync only
            // the reduce stream (near-zero cost if physics integration has run since launch)
            if(!m_hasAsyncResult)
                return true;  // Conservative on first iteration: no result available yet
            grains::cuda::synchronizeLinkedCellReduction(*m_reductionBackend.state);
            return (T(4) * m_maxDisplacementSquared[0] > m_skinThickness * m_skinThickness);
        }
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Generates neighbor cells. */
    void generateNeighborCells()
    {
        GAssert(m_numCells > 0, "LinkedCell::generateNeighborCells: number of cells is zero!");

        if constexpr(M == MemType::HOST)
        {
            m_cells[0]->generateNeighborCells(m_neighborCells.getData());
        }
        else if constexpr(M == MemType::DEVICE)
        {
            grains::cuda::generateNeighborCells(m_cells.getData(),
                                                m_numCells,
                                                m_neighborCells.getData());
        }
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Links obstacles to cells. */
    void linkObstacles()
    {
        // If there is no obstacle
        if(m_numObstacles == 0)
            return;

        if constexpr(M == MemType::HOST)
        {
            // Lambda to extract support point from rigid body in given world direction
            auto support
                = [this](uint obstacleIndex, const Vector3<T>& worldDirection) -> Vector3<T> {
                // Transform world direction to local coordinates using inverse rotation
                const Quaternion<T>& q              = m_quaternions->at(obstacleIndex);
                const Vector3<T>     localDirection = q << worldDirection;
                Vector3<T> supPt = (*m_rb)[obstacleIndex]->getConvex()->support(localDirection);
                transform(q, m_positions->at(obstacleIndex), supPt);
                return supPt;
            };

            // Cell info
            const CellDimensions numCells = m_cells[0]->getNumCellsPerDirection();

            for(uint i = 0; i < m_numObstacles; ++i)
            {
                // offset in the obstacleCellID buffer
                const uint offset = i * m_maxCellsPerObstacle;

                // AABB by querying support in all 6 axis directions
                const Vector3<T> minExt(support(i, Vector3<T>(-1, 0, 0))[X],
                                        support(i, Vector3<T>(0, -1, 0))[Y],
                                        support(i, Vector3<T>(0, 0, -1))[Z]);
                const Vector3<T> maxExt(support(i, Vector3<T>(1, 0, 0))[X],
                                        support(i, Vector3<T>(0, 1, 0))[Y],
                                        support(i, Vector3<T>(0, 0, 1))[Z]);

                // Convert world coordinates to cell coordinates
                const CellIndex minCell = m_cells[0]->computeCellID(minExt, false);
                const CellIndex maxCell = m_cells[0]->computeCellID(maxExt, false);

                uint cellCount = 0;
                int  minX      = std::max((int)minCell.x - 1, 0);
                int  maxX      = std::min((int)maxCell.x + 1, (int)numCells.x - 1);
                int  minY      = std::max((int)minCell.y - 1, 0);
                int  maxY      = std::min((int)maxCell.y + 1, (int)numCells.y - 1);
                int  minZ      = std::max((int)minCell.z - 1, 0);
                int  maxZ      = std::min((int)maxCell.z + 1, (int)numCells.z - 1);

                // Nested loops with 1-ring expansion; guard against overflow
                for(int x = minX; x <= maxX && cellCount < m_maxCellsPerObstacle; ++x)
                {
                    for(int y = minY; y <= maxY && cellCount < m_maxCellsPerObstacle; ++y)
                    {
                        for(int z = minZ; z <= maxZ && cellCount < m_maxCellsPerObstacle; ++z)
                        {
                            uint cellHash = m_cells[0]->computeCellHash(
                                makeCellIndex((uint)x, (uint)y, (uint)z));
                            m_obstacleCellID[offset + cellCount] = cellHash;
                            ++cellCount;
                        }
                    }
                }

                // Update the count in obstacleID buffer
                m_obstacleID[i].x = i;
                m_obstacleID[i].y = cellCount;
            }
        }
        else if constexpr(M == MemType::DEVICE)
        {
            grains::cuda::linkObstacles(m_rb->getData(),
                                        m_positions->getData(),
                                        m_quaternions->getData(),
                                        m_cells.getData(),
                                        m_numObstacles,
                                        m_maxCellsPerObstacle,
                                        m_obstacleID.getData(),
                                        m_obstacleCellID.getData());
        }

        // Obstacles are now in sync; reset flag so we don't relink every iteration for
        // static obstacles. Obstacle time integration must set this back to true when
        // obstacles actually move.
        GrainsParameters<T>::m_simulationState.obstaclesMoved = false;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Computes the maximum displacement. */
    T computeMaxDisplacement()
    {
        if constexpr(M == MemType::HOST)
        {
            T maxDisplacementSquared = 0;
            for(uint i = 0; i < m_oldPosition.getSize(); ++i)
            {
                T dispSquared = norm2(m_positions->at(i) - m_oldPosition[i]);
                if(dispSquared > maxDisplacementSquared)
                    maxDisplacementSquared = dispSquared;
            }
            m_maxDisplacementSquared[0] = maxDisplacementSquared;
        }
        else if constexpr(M == MemType::DEVICE)
        {
            grains::cuda::computeMaxDisplacementSquared(*m_reductionBackend.state,
                                                        m_oldPosition.getData(),
                                                        m_positions->getData(),
                                                        m_positions->getSize(),
                                                        m_maxDisplacementSquared.getData());
        }

        return m_maxDisplacementSquared[0];
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Computes the skin thickness based on the maximum displacement. */
    // ---------------------------------------------------------------------------------------------
    /** @brief Launches the displacement reduce asynchronously on the dedicated reduce stream.
        Called at the end of each updateLinkedCells() so the result is ready for the next
        iteration's needsUpdate() check with near-zero synchronization cost. */
    void launchDisplacementReduceAsync()
    {
        if constexpr(M == MemType::DEVICE)
        {
            grains::cuda::launchMaxDisplacementSquared(*m_reductionBackend.state,
                                                       m_oldPosition.getData(),
                                                       m_positions->getData(),
                                                       m_positions->getSize(),
                                                       m_maxDisplacementSquared.getData());
            m_hasAsyncResult = true;
        }
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Computes the skin thickness based on the maximum displacement. */
    T computeSkinThickness() const
    {
        // Smoothing factor
        constexpr T mu = T(0.4);
        // Max Cap the skin thickness at 20% of the cell size
        constexpr T maxSkinThickness = T(0.2);

        // Guard: if called before any iteration has elapsed (e.g. sort fires on iteration 0),
        // keep the current skin to avoid division by zero producing NaN.
        if(m_numIterationsSinceLastUpdate == 0)
            return m_skinThickness;

        const T newThickness = T(2) * sqrt(m_maxDisplacementSquared[0]) * m_updateFrequency
                               / m_numIterationsSinceLastUpdate;
        T skinThickness = mu * newThickness + (1 - mu) * m_skinThickness;
        if(skinThickness > maxSkinThickness * m_cellSizeWithoutSkin)
            skinThickness = maxSkinThickness * m_cellSizeWithoutSkin;

        return skinThickness;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Updates the linked cells and returns whether the LinkedCell has been updated. */
    virtual bool updateLinkedCells() = 0;
    //@}
};

#endif
