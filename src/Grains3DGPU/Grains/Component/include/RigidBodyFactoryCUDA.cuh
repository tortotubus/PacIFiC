#ifndef RIGID_BODY_FACTORY_CUDA_CUH
#define RIGID_BODY_FACTORY_CUDA_CUH

#include "GrainsMemBuffer.hh"
#include "RigidBody.hh"

namespace grains::cuda
{
    /** Construct device-side rigid bodies from host-side polymorphic objects. */
    template <typename T>
    void copyRigidBodiesToDevice(GrainsMemBuffer<RigidBody<T>*, MemType::HOST>&   hostBodies,
                                 GrainsMemBuffer<RigidBody<T>*, MemType::DEVICE>& deviceBodies);

    /** Delete rigid bodies previously constructed on the device. */
    template <typename T>
    void freeRigidBodiesOnDevice(
        GrainsMemBuffer<RigidBody<T>*, MemType::DEVICE>& deviceBodies);
}

#endif
