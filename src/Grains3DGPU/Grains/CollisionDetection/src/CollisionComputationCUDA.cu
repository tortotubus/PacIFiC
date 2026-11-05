#include "CollisionComputationCUDA.hh"
#include "CollisionComputation.hh"
#include "CollisionDispatch.hh"
#include "GrainsCUDAUtils.hh"
#include "GrainsParameters.hh"

namespace
{
    template <typename T>
    __global__ void relativeKernel(RelativeTransformView<T> view)
    {
        view.compute(blockIdx.x * blockDim.x + threadIdx.x);
    }
    template <GJKType Variant, bool Accelerated, BoundingVolumeType Bounds, typename View>
    __global__ void detectKernel(View view, bool prebuilt)
    {
        view.template detect<Variant, Accelerated, Bounds>(blockIdx.x * blockDim.x + threadIdx.x,
                                                           prebuilt);
    }
    template <typename T>
    __global__ void transformKernel(ContactTransformView<T> view)
    {
        view.transform(blockIdx.x * blockDim.x + threadIdx.x);
    }
    template <typename T>
    void relative(RelativeTransformView<T> view)
    {
        if(view.count == 0)
            return;
        uint blocks, threads;
        computeOptimalThreadsAndBlocks(view.count, GrainsParameters<T>::m_GPU, blocks, threads);
        relativeKernel<<<blocks, threads>>>(view);
        cudaErrCheck(cudaGetLastError());
        cudaErrCheck(cudaDeviceSynchronize());
    }
    template <typename T, typename View>
    void detect(View view, NarrowPhaseOptions options)
    {
        if(view.data.count == 0)
            return;
        uint blocks, threads;
        computeOptimalThreadsAndBlocks(view.data.count,
                                       GrainsParameters<T>::m_GPU,
                                       blocks,
                                       threads);
        grains::detail::dispatchCollision(options,
                                          [&](auto variant, auto accelerated, auto bounds) {
                                              detectKernel<decltype(variant)::value,
                                                           decltype(accelerated)::value,
                                                           decltype(bounds)::value>
                                                  <<<blocks, threads>>>(view,
                                                                        options.prebuiltShapes);
                                          });
        cudaErrCheck(cudaGetLastError());
        cudaErrCheck(cudaDeviceSynchronize());
    }
    template <typename T>
    void transform(ContactTransformView<T> view)
    {
        if(view.count == 0)
            return;
        uint blocks, threads;
        computeOptimalThreadsAndBlocks(view.count, GrainsParameters<T>::m_GPU, blocks, threads);
        transformKernel<<<blocks, threads>>>(view);
        cudaErrCheck(cudaGetLastError());
        cudaErrCheck(cudaDeviceSynchronize());
    }
}
namespace grains::cuda
{
    void computeRelativeTransforms(RelativeTransformView<float> view)
    {
        relative(view);
    }
    void computeRelativeTransforms(RelativeTransformView<double> view)
    {
        relative(view);
    }
    void detectRelativeCollisions(RelativeCollisionView<float> view, NarrowPhaseOptions options)
    {
        detect<float>(view, options);
    }
    void detectRelativeCollisions(RelativeCollisionView<double> view, NarrowPhaseOptions options)
    {
        detect<double>(view, options);
    }
    void detectGlobalCollisions(GlobalCollisionView<float> view, NarrowPhaseOptions options)
    {
        detect<float>(view, options);
    }
    void detectGlobalCollisions(GlobalCollisionView<double> view, NarrowPhaseOptions options)
    {
        detect<double>(view, options);
    }
    void transformContacts(ContactTransformView<float> view)
    {
        transform(view);
    }
    void transformContacts(ContactTransformView<double> view)
    {
        transform(view);
    }
}
