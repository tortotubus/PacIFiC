#ifndef _LINKEDCELLFACTORY_HH_
#define _LINKEDCELLFACTORY_HH_

#include "GrainsMemBuffer.hh"
#include "GrainsParameters.hh"
#include "LinkedCell.hh"
#include "LinkedCellFactoryCUDA.cuh"
#include "LinkedCell_Host.hh"

// =================================================================================================
/** @brief The class LinkedCellFactory.

    Creates the linked cell structure for the simulation. Host construction stays in this ordinary
    C++ header; device implementation selection is owned by the CUDA backend.

    @author A.YAZDANI - 2025 - Construction */
// =================================================================================================
template <typename T, MemType M>
class LinkedCellFactory
{
    static_assert(M == MemType::HOST || M == MemType::DEVICE,
                  "LinkedCellFactory only supports MemType::HOST or MemType::DEVICE");

private:
    /** @name Constructors & Destructor */
    //@{
    /** @brief Default constructor (forbidden) */
    LinkedCellFactory() = default;

    /** @brief Destructor (forbidden) */
    ~LinkedCellFactory() = default;
    //@}

public:
    /**@name Methods */
    //@{
    /** @brief Creates the configured linked-cell implementation.
        @param rigidBodies Rigid-body buffer
        @param positions Position buffer
        @param orientations Quaternion buffer
        @param parameters Linked-cell parameters
        @param numObstacles Number of obstacles
        @param numParticles Number of particles */
    static std::unique_ptr<LinkedCell<T, M>>
        create(const GrainsMemBuffer<RigidBody<T>*, M>* rigidBodies,
               const GrainsMemBuffer<Vector3<T>, M>&    positions,
               const GrainsMemBuffer<Quaternion<T>, M>& orientations,
               const LinkedCellParameters<T>&           parameters,
               uint                                     numObstacles,
               uint                                     numParticles)
    {
        if constexpr(M == MemType::HOST)
        {
            return std::make_unique<LinkedCell_Host<T>>(rigidBodies,
                                                        positions,
                                                        orientations,
                                                        parameters,
                                                        numObstacles,
                                                        numParticles);
        }
        else
        {
            return grains::cuda::createLinkedCell<T>(rigidBodies,
                                                     positions,
                                                     orientations,
                                                     parameters,
                                                     numObstacles,
                                                     numParticles);
        }
    }
    //@}
};

#endif
