#include "ContactForceModelFactoryCUDA.cuh"

#include "GrainsCUDAError.hh"
#include "GrainsParameters.hh"
#include "HookeContactForceModel.hh"
#include "HookeMemoryContactForceModel.hh"

namespace grains::cuda
{
    template <typename T, typename... Arguments>
    __GLOBAL__ void
        createContactForceModelKernel(ContactForceModel<T>** models, uint index, Arguments... args)
    {
        const uint threadId = blockIdx.x * blockDim.x + threadIdx.x;
        if(threadId > 0)
            return;

        if constexpr(sizeof...(args) == 5)
        {
            models[index] = new HookeContactForceModel<T>(args...);
        }
        else if constexpr(sizeof...(args) == 8)
        {
            models[index] = new HookeMemoryContactForceModel<T>(args...);
        }
        else
        {
            static_assert(sizeof...(args) == 5 || sizeof...(args) == 8,
                          "Unsupported contact-model parameter count");
        }
    }

    template <typename T>
    void copyContactForceModelsToDevice(
        GrainsMemBuffer<ContactForceModel<T>*, MemType::HOST>&   hostModels,
        GrainsMemBuffer<ContactForceModel<T>*, MemType::DEVICE>& deviceModels)
    {
        deviceModels.initialize(hostModels.getSize());
        for(uint i = 0; i < hostModels.getSize(); ++i)
        {
            if(hostModels[i] == nullptr)
                continue;

            ContactForceModel<T>* contact     = hostModels[i];
            ContactForceModelType contactType = contact->getContactForceModelType();

            if(contactType == HOOKE)
            {
                auto* model = dynamic_cast<HookeContactForceModel<T>*>(contact);
                T     kn, en, etat, muc, kr;
                model->getContactForceModelParameters(kn, en, etat, muc, kr);
                createContactForceModelKernel<<<1, 1>>>(deviceModels.getData(),
                                                        i,
                                                        kn,
                                                        en,
                                                        etat,
                                                        muc,
                                                        kr);
            }
            else if(contactType == HOOKEMEMORY)
            {
                auto* model = dynamic_cast<HookeMemoryContactForceModel<T>*>(contact);
                T     kn, en, kt, etat, muc, mur, etarpf;
                model->getContactForceModelParameters(kn, en, kt, etat, muc, mur, etarpf);
                createContactForceModelKernel<<<1, 1>>>(deviceModels.getData(),
                                                        i,
                                                        kn,
                                                        en,
                                                        kt,
                                                        etat,
                                                        muc,
                                                        mur,
                                                        etarpf,
                                                        GrainsParameters<T>::m_dt);
            }
            else
            {
                GAbort("Contact force model is not implemented for GPU!", "Aborting Grains!");
            }
        }
        cudaErrCheck(cudaDeviceSynchronize());
    }

    template void copyContactForceModelsToDevice(
        GrainsMemBuffer<ContactForceModel<float>*, MemType::HOST>&,
        GrainsMemBuffer<ContactForceModel<float>*, MemType::DEVICE>&);
    template void copyContactForceModelsToDevice(
        GrainsMemBuffer<ContactForceModel<double>*, MemType::HOST>&,
        GrainsMemBuffer<ContactForceModel<double>*, MemType::DEVICE>&);
}  // namespace grains::cuda
