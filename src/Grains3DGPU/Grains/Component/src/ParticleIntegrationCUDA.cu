#include "ParticleIntegrationCUDA.hh"
#include "ParticleIntegration.hh"
#include "GrainsCUDAUtils.hh"
#include "GrainsParameters.hh"

namespace
{
    template <typename T>
    __global__ void moveKernel(ParticleIntegrationView<T> view)
    {
        const uint offset = blockIdx.x * blockDim.x + threadIdx.x;
        if(offset < view.end - view.begin)
            view.move(view.begin + offset);
    }

    template <typename T>
    __global__ void advanceKernel(ParticleVelocityView<T> view)
    {
        const uint offset = blockIdx.x * blockDim.x + threadIdx.x;
        if(offset < view.end - view.begin)
            view.advanceVelocity(view.begin + offset);
    }

    template <typename T>
    __global__ void updateKernel(SubBodyTransformView<T> view)
    {
        const uint offset = blockIdx.x * blockDim.x + threadIdx.x;
        if(offset < view.end - view.begin)
            view.update(view.begin + offset, EPS<T>);
    }

    template <typename T>
    void move(ParticleIntegrationView<T> view)
    {
        if(view.begin >= view.end)
            return;
        uint blocks, threads;
        computeOptimalThreadsAndBlocks(view.end - view.begin,
                                       GrainsParameters<T>::m_GPU,
                                       blocks,
                                       threads);
        moveKernel<<<blocks, threads>>>(view);
        cudaErrCheck(cudaGetLastError());
    }

    template <typename T>
    void advance(ParticleVelocityView<T> view)
    {
        if(view.begin >= view.end)
            return;
        uint blocks, threads;
        computeOptimalThreadsAndBlocks(view.end - view.begin,
                                       GrainsParameters<T>::m_GPU,
                                       blocks,
                                       threads);
        advanceKernel<<<blocks, threads>>>(view);
        cudaErrCheck(cudaGetLastError());
    }

    template <typename T>
    void update(SubBodyTransformView<T> view)
    {
        if(view.begin >= view.end)
            return;
        uint blocks, threads;
        computeOptimalThreadsAndBlocks(view.end - view.begin,
                                       GrainsParameters<T>::m_GPU,
                                       blocks,
                                       threads);
        updateKernel<<<blocks, threads>>>(view);
        cudaErrCheck(cudaGetLastError());
    }
}

namespace grains::cuda
{
    void moveParticles(ParticleIntegrationView<float> view)
    {
        move(view);
    }
    void moveParticles(ParticleIntegrationView<double> view)
    {
        move(view);
    }
    void advanceVelocity(ParticleVelocityView<float> view)
    {
        advance(view);
    }
    void advanceVelocity(ParticleVelocityView<double> view)
    {
        advance(view);
    }
    void updateSubBodyPositions(SubBodyTransformView<float> view)
    {
        update(view);
    }
    void updateSubBodyPositions(SubBodyTransformView<double> view)
    {
        update(view);
    }
}
