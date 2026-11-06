#include "TimeIntegratorFactoryCUDA.cuh"

#include "FirstOrderExplicit.hh"
#include "GrainsCUDAError.hh"
#include "SecondOrderLeapFrog.hh"

namespace grains::cuda
{
    template <typename T>
    __GLOBAL__ void createTimeIntegratorKernel(TimeIntegrator<T>**      integrators,
                                               const uint               index,
                                               const TimeIntegratorType type,
                                               const T                  dt)
    {
        const uint threadId = blockIdx.x * blockDim.x + threadIdx.x;
        if(threadId > 0)
            return;

        if(type == FIRSTORDEREXPLICIT)
            integrators[index] = new FirstOrderExplicit<T>(dt);
        else if(type == SECONDORDERLEAPFROG)
            integrators[index] = new SecondOrderLeapFrog<T>(dt);
        else
            GAbort("Time integrator is not implemented for GPU!", "Aborting Grains!");
    }

    template <typename T>
    void copyTimeIntegratorsToDevice(
        GrainsMemBuffer<TimeIntegrator<T>*, MemType::HOST>&   hostIntegrators,
        GrainsMemBuffer<TimeIntegrator<T>*, MemType::DEVICE>& deviceIntegrators)
    {
        deviceIntegrators.initialize(hostIntegrators.getSize());
        for(uint i = 0; i < hostIntegrators.getSize(); ++i)
        {
            if(hostIntegrators[i] == nullptr)
                continue;

            createTimeIntegratorKernel<<<1, 1>>>(deviceIntegrators.getData(),
                                                 i,
                                                 hostIntegrators[i]->getTimeIntegratorType(),
                                                 hostIntegrators[i]->getTimeStep());
        }
        cudaErrCheck(cudaDeviceSynchronize());
    }

    template void
        copyTimeIntegratorsToDevice(GrainsMemBuffer<TimeIntegrator<float>*, MemType::HOST>&,
                                    GrainsMemBuffer<TimeIntegrator<float>*, MemType::DEVICE>&);
    template void
        copyTimeIntegratorsToDevice(GrainsMemBuffer<TimeIntegrator<double>*, MemType::HOST>&,
                                    GrainsMemBuffer<TimeIntegrator<double>*, MemType::DEVICE>&);
}  // namespace grains::cuda
