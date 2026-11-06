#include "LinkedCellFactoryCUDA.cuh"

#include "GrainsAssert.hh"
#include "GrainsParameters.hh"
#include "LinkedCell_Atomic.cuh"
#include "LinkedCell_AtomicFixed.cuh"
#include "LinkedCell_SortBased.cuh"

namespace grains::cuda
{
    template <typename T>
    std::unique_ptr<LinkedCell<T, MemType::DEVICE>>
        createLinkedCell(const GrainsMemBuffer<RigidBody<T>*, MemType::DEVICE>* rigidBodies,
                         const GrainsMemBuffer<Vector3<T>, MemType::DEVICE>&    positions,
                         const GrainsMemBuffer<Quaternion<T>, MemType::DEVICE>& orientations,
                         const LinkedCellParameters<T>&                         parameters,
                         unsigned int                                           numObstacles,
                         unsigned int                                           numParticles)
    {
        switch(parameters.type)
        {
        case LinkedCellType::SORTBASED:
            return std::make_unique<LinkedCell_SortBased<T>>(rigidBodies,
                                                             positions,
                                                             orientations,
                                                             parameters,
                                                             numObstacles,
                                                             numParticles);
        case LinkedCellType::ATOMIC:
            return std::make_unique<LinkedCell_Atomic<T>>(rigidBodies,
                                                          positions,
                                                          orientations,
                                                          parameters,
                                                          numObstacles,
                                                          numParticles);
        case LinkedCellType::ATOMICFIXED:
            return std::make_unique<LinkedCell_AtomicFixed<T>>(rigidBodies,
                                                               positions,
                                                               orientations,
                                                               parameters,
                                                               numObstacles,
                                                               numParticles);
        case LinkedCellType::HOST:
            GAbort("Host linked cells cannot be constructed by the CUDA backend!");
        }

        return nullptr;
    }

    template std::unique_ptr<LinkedCell<float, MemType::DEVICE>>
        createLinkedCell(const GrainsMemBuffer<RigidBody<float>*, MemType::DEVICE>*,
                         const GrainsMemBuffer<Vector3<float>, MemType::DEVICE>&,
                         const GrainsMemBuffer<Quaternion<float>, MemType::DEVICE>&,
                         const LinkedCellParameters<float>&,
                         unsigned int,
                         unsigned int);
    template std::unique_ptr<LinkedCell<double, MemType::DEVICE>>
        createLinkedCell(const GrainsMemBuffer<RigidBody<double>*, MemType::DEVICE>*,
                         const GrainsMemBuffer<Vector3<double>, MemType::DEVICE>&,
                         const GrainsMemBuffer<Quaternion<double>, MemType::DEVICE>&,
                         const LinkedCellParameters<double>&,
                         unsigned int,
                         unsigned int);
}  // namespace grains::cuda
