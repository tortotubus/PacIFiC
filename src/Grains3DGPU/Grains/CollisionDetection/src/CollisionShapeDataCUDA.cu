#include "CollisionShapeData.hh"
#include "CollisionShapeDataCUDA.hh"
#include "GrainsCUDAUtils.hh"
#include "GrainsMemBuffer.hh"
#include "GrainsParameters.hh"

namespace
{
    template <typename T, bool Bounds>
    __global__ void fillKernel(ShapeTableView<T> view, const uint* representatives)
    {
        const uint i = blockIdx.x * blockDim.x + threadIdx.x;
        if constexpr(Bounds)
            view.fillBounds(i, representatives);
        else
            view.fillShape(i, representatives);
    }
    template <typename T>
    void build(ShapeTableView<T> view, const uint* representatives)
    {
        if(view.count == 0)
            return;
        GrainsMemBuffer<uint, MemType::DEVICE> deviceRepresentatives;
        deviceRepresentatives.initialize(view.count);
        cudaErrCheck(cudaMemcpy(deviceRepresentatives.getData(),
                                representatives,
                                view.count * sizeof(uint),
                                cudaMemcpyHostToDevice));
        uint blocks, threads;
        computeOptimalThreadsAndBlocks(view.count, GrainsParameters<T>::m_GPU, blocks, threads);
        fillKernel<T, false><<<blocks, threads>>>(view, deviceRepresentatives.getData());
        cudaErrCheck(cudaGetLastError());
        fillKernel<T, true><<<blocks, threads>>>(view, deviceRepresentatives.getData());
        cudaErrCheck(cudaGetLastError());
        cudaErrCheck(cudaDeviceSynchronize());
    }
}
namespace grains::cuda
{
    void buildShapeTables(ShapeTableView<float> view, const uint* representatives)
    {
        build(view, representatives);
    }
    void buildShapeTables(ShapeTableView<double> view, const uint* representatives)
    {
        build(view, representatives);
    }
}
