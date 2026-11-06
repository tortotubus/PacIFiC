#ifndef GRAINS_LINKED_CELL_SORT_BASED_CUH
#define GRAINS_LINKED_CELL_SORT_BASED_CUH

#include "GrainsPlatform.hh"

#include "GrainsMemBuffer.hh"
#include "LinkedCell.hh"

// =================================================================================================
/** @brief The class LinkedCell_SortBased.

    This class provides functionalities to manage linked cells for collision detection in the
    simulation using a sort-based approach. It is a derived class of LinkedCell and implements the
    update of linked cells based on sorting the particle hashes. This is designed to work on the
    device (GPU). This has optimal space complexity, while the time complexity is O(n) for computing
    the particle hashes and O(nk) for sorting the particle IDs using a parallel radix sort
    algorithm, where n is the number of particles and k is the number of cells.

    @author A.Yazdani - 2025 - Construction */
// =================================================================================================
template <typename T>
class LinkedCell_SortBased : public LinkedCell<T, MemType::DEVICE>
{
    using LC = LinkedCell<T, MemType::DEVICE>;
    using LC::m_cells;
    using LC::m_cellSizeWithoutSkin;
    using LC::m_maxDisplacementSquared;
    using LC::m_neighborCells;
    using LC::m_numCells;
    using LC::m_numIterationsSinceLastUpdate;
    using LC::m_numObstacles;
    using LC::m_numParticles;
    using LC::m_oldPosition;
    using LC::m_positions;
    using LC::m_skinThickness;
    using LC::m_useAdaptiveSkin;

protected:
    /** @name Parameters */
    //@{
    /** \brief CUDA-owned streams, event, workspace, and SortBased output buffers. */
    grains::cuda::SortBasedStatePtr m_cudaState;
    //@}

public:
    /** @name Constructors */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Default constructor */
    LinkedCell_SortBased() = default;

    // ---------------------------------------------------------------------------------------------
    /** @brief Constructor with parameters
        @param rb Rigid body buffer
        @param positions Positions buffer
        @param quaternions Quaternions buffer
        @param linkedCellParameters Linked cell parameters
        @param nObstacles number of obstacles
        @param nParticles number of particles */
    LinkedCell_SortBased(const GrainsMemBuffer<RigidBody<T>*, MemType::DEVICE>* rb,
                         const GrainsMemBuffer<Vector3<T>, MemType::DEVICE>&    positions,
                         const GrainsMemBuffer<Quaternion<T>, MemType::DEVICE>& quaternions,
                         const LinkedCellParameters<T>& linkedCellParameters,
                         const uint                     nObstacles,
                         const uint                     nParticles)
        : LinkedCell<T, MemType::DEVICE>(
              rb, positions, quaternions, linkedCellParameters, nObstacles, nParticles)
        , m_cudaState(grains::cuda::createSortBasedState(nParticles, m_numCells))
    {
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Destructor */
    ~LinkedCell_SortBased() override = default;
    //@}

    /** @name Get methods */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Gets packed cell-particle IDs. */
    const uint64_t* getCellParticleIDs() const override
    {
        return m_cudaState ? grains::cuda::sortBasedCellParticleIDs(*m_cudaState) : nullptr;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Gets cell prefix sums (start indices for each cell) */
    const uint* getCellPrefixSums() const override
    {
        return m_cudaState ? grains::cuda::sortBasedCellPrefixSums(*m_cudaState) : nullptr;
    }
    //@}

    /** @name Methods */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Prepares packed particle IDs and optionally resizes the adaptive grid. */
    void prepareLinkedCellUpdate()
    {
        T cellSize = m_cellSizeWithoutSkin;
        if(m_useAdaptiveSkin)
        {
            m_skinThickness                = this->computeSkinThickness();
            cellSize                       = m_cellSizeWithoutSkin + m_skinThickness;
            m_maxDisplacementSquared[0]    = T(0);
            m_numIterationsSinceLastUpdate = 0;
        }

        grains::cuda::SortBasedUpdateView<T> view;
        view.cells         = m_cells.getData();
        view.positions     = m_positions->getData();
        view.oldPositions  = m_oldPosition.getData();
        view.neighborCells = &m_neighborCells;
        view.numCells      = m_numCells;
        view.numObstacles  = m_numObstacles;
        view.numParticles  = m_numParticles;
        view.adaptive      = m_useAdaptiveSkin;
        view.cellSize      = cellSize;
        m_numCells         = grains::cuda::prepareSortBasedUpdate(*m_cudaState, view);
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Updates the linked cells based on the transformations */
    bool updateLinkedCells() override
    {
        auto& SS = GrainsParameters<T>::m_simulationState;

        if(this->needsUpdate())
        {
            prepareLinkedCellUpdate();
            // Relink obstacles
            this->linkObstacles();
        }
        else
        {
            // Relink obstacles if they moved
            if(SS.obstaclesMoved)
                this->linkObstacles();
            if(m_useAdaptiveSkin)
                this->launchDisplacementReduceAsync();
            return false;
        }

        grains::cuda::finalizeSortBasedUpdate<T>(*m_cudaState, m_numParticles);

        // Launch displacement reduce async so next iteration's check costs almost nothing
        if(m_useAdaptiveSkin)
            this->launchDisplacementReduceAsync();

        return true;
    }
};

#endif
