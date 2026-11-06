#ifndef TIME_INTEGRATOR_FACTORY_CUDA_CUH
#define TIME_INTEGRATOR_FACTORY_CUDA_CUH

#include "GrainsMemBuffer.hh"
#include "TimeIntegrator.hh"

namespace grains::cuda
{
    /** Construct device-side integrators from the corresponding host-side polymorphic objects. */
    template <typename T>
    void copyTimeIntegratorsToDevice(
        GrainsMemBuffer<TimeIntegrator<T>*, MemType::HOST>&   hostIntegrators,
        GrainsMemBuffer<TimeIntegrator<T>*, MemType::DEVICE>& deviceIntegrators);
}

#endif
