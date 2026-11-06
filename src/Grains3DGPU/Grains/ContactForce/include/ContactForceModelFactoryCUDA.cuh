#ifndef CONTACT_FORCE_MODEL_FACTORY_CUDA_CUH
#define CONTACT_FORCE_MODEL_FACTORY_CUDA_CUH

#include "ContactForceModel.hh"
#include "GrainsMemBuffer.hh"

namespace grains::cuda
{
    /** Construct device-side contact models from the corresponding host-side polymorphic objects. */
    template <typename T>
    void copyContactForceModelsToDevice(
        GrainsMemBuffer<ContactForceModel<T>*, MemType::HOST>&   hostModels,
        GrainsMemBuffer<ContactForceModel<T>*, MemType::DEVICE>& deviceModels);
}

#endif
