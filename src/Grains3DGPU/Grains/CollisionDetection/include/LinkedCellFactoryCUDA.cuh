#ifndef LINKED_CELL_FACTORY_CUDA_CUH
#define LINKED_CELL_FACTORY_CUDA_CUH

#include <memory>

#include "GrainsMemBuffer.hh"

template <typename T>
class Quaternion;
template <typename T>
class RigidBody;
template <typename T>
class Vector3;
template <typename T, MemType M>
class LinkedCell;
template <typename T>
struct LinkedCellParameters;

namespace grains::cuda
{
    /** Construct the configured host-side manager for a device linked-cell implementation. */
    template <typename T>
    std::unique_ptr<LinkedCell<T, MemType::DEVICE>>
        createLinkedCell(const GrainsMemBuffer<RigidBody<T>*, MemType::DEVICE>* rigidBodies,
                         const GrainsMemBuffer<Vector3<T>, MemType::DEVICE>&    positions,
                         const GrainsMemBuffer<Quaternion<T>, MemType::DEVICE>& orientations,
                         const LinkedCellParameters<T>&                         parameters,
                         unsigned int                                           numObstacles,
                         unsigned int                                           numParticles);
}

#endif
