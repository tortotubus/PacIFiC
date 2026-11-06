#ifndef _NEIGHBORLIST_LINKEDCELL_HH_
#define _NEIGHBORLIST_LINKEDCELL_HH_

#include "GrainsPlatform.hh"

#include "GrainsMemBuffer.hh"
#include "GrainsParameters.hh"
#include "LinkedCell.hh"
#include "LinkedCellFactory.hh"
#include "LinkedCell_Host.hh"
#include "NeighborList.hh"
#include "NeighborListLinkedCellCUDA.hh"
#include "NeighborList_Kernels.hh"

namespace grains::detail
{
    template <MemType M>
    struct NeighborListLinkedCellStorage
    {
    };

    template <>
    struct NeighborListLinkedCellStorage<MemType::DEVICE>
    {
        grains::cuda::NeighborListLinkedCellStatePtr state;
    };
}

// =================================================================================================
/** @brief The class NeighborList_LinkedCell.

    This is a derived class of NeighborList. It implements the neighbor list creation using an O(n)
    algorithm. This is useful for systems with large number of components.

    @author A.Yazdani - 2025 - Construction */
// =================================================================================================
template <typename T, MemType M>
class NeighborList_LinkedCell : public NeighborList<T, M>
{
    using NL = NeighborList<T, M>;
    using NL::m_pairCount;
    using NL::m_pairList;

protected:
    /** @name Parameters */
    //@{
    /** \brief LinkedCell */
    std::unique_ptr<LinkedCell<T, M>> m_LinkedCell;
    /** \brief CUDA-owned buffers, scan workspace, counter, and streams. */
    grains::detail::NeighborListLinkedCellStorage<M> m_cudaBackend;
    //@}

public:
    /** @name Constructors */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Constructor */
    NeighborList_LinkedCell() = default;

    // ---------------------------------------------------------------------------------------------
    /** @brief Constructor with parameters
        @param rb Rigid body buffer
        @param positions Positions buffer
        @param quaternions Quaternions buffer
        @param linkedCellParameters Linked cell parameters
        @param nObstacles number of obstacles
        @param nParticles number of particles */
    NeighborList_LinkedCell(const GrainsMemBuffer<RigidBody<T>*, M>* rb,
                            const GrainsMemBuffer<Vector3<T>, M>&    positions,
                            const GrainsMemBuffer<Quaternion<T>, M>& quaternions,
                            const LinkedCellParameters<T>&           linkedCellParameters,
                            const uint                               nObstacles,
                            const uint                               nParticles)
    {
        // Create the LinkedCell buffer
        m_LinkedCell = LinkedCellFactory<T, M>::create(rb,
                                                       positions,
                                                       quaternions,
                                                       linkedCellParameters,
                                                       nObstacles,
                                                       nParticles);

        // Initialize pair list
        auto initNumPairs = GrainsParameters<T>::m_collisionDetection.linkedCellParameters
                                .initialNumberOfPairsPerParticle;
        m_pairList.initialize(initNumPairs * nParticles);
        m_pairList.fill();

        if constexpr(M == MemType::DEVICE)
        {
            m_cudaBackend.state = grains::cuda::createNeighborListLinkedCellState(nParticles);
        }
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Destructor */
    ~NeighborList_LinkedCell() override = default;
    //@}

    /** @name Methods */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Updates the neighbor list
        @param positions array of positions
        @param nObstacles number of obstacles
        @param nParticles number of particles */
    bool updateNeighborList(GrainsMemBuffer<Vector3<T>, M>& positions,
                            const uint                      nObstacles,
                            const uint                      nParticles) final
    {
        // Update linked cells (internally uses streams for GPU variants)
        bool LC_updated = m_LinkedCell->updateLinkedCells();

        if(LC_updated)
        {
            if constexpr(M == MemType::HOST)
            {
                auto* LC_host = static_cast<LinkedCell_Host<T>*>(m_LinkedCell.get());
                updateNeighborList_LC_Host(LC_host->getCellNeighborsList(),
                                           LC_host->getObstacleIDs(),
                                           LC_host->getObstacleCellIDs(),
                                           LC_host->getParticleIDs(),
                                           LC_host->getCellIDs(),
                                           LC_host->getCellParticles(),
                                           LC_host->getMaxCellsPerObstacle(),
                                           nObstacles,
                                           nParticles,
                                           m_pairList);
                m_pairCount[0] = m_pairList.getSize();
            }
            else if constexpr(M == MemType::DEVICE)
            {
                grains::cuda::NeighborListLinkedCellUpdateView<T> view;
                view.linkedCell   = m_LinkedCell.get();
                view.pairList     = &m_pairList;
                view.pairCount    = m_pairCount.getData();
                view.numObstacles = nObstacles;
                view.numParticles = nParticles;
                grains::cuda::updateLinkedCellNeighborList(*m_cudaBackend.state, view);
            }

            return true;
        }
        else
            return false;
    }
    //@}
};

#endif
