#include "CellsFactory.hh"
#include "GrainsCUDAUtils.hh"
#include "ParticleSortCUDA.hh"
#include "ParticleSortingCUDA.hh"
#include "ParticleSortingStorage.hh"
#include <array>

namespace
{
    template <typename T>
    __global__ void deleteCells(Cells<T, CellOrdering::MORTON>** cells)
    {
        delete cells[0];
        cells[0] = nullptr;
    }
    template <typename T>
    __global__ void mortonKernel(ParticleMortonView<T>                  view,
                                 Cells<T, CellOrdering::MORTON>* const* cells)
    {
        view.compute(blockIdx.x * blockDim.x + threadIdx.x, **cells);
    }
    template <typename Value>
    __global__ void gatherKernel(IndexedGatherView<Value> view)
    {
        view.gather(blockIdx.x * blockDim.x + threadIdx.x);
    }
    __global__ void masterSlotKernel(MasterSlotView view)
    {
        view.rebuild(blockIdx.x * blockDim.x + threadIdx.x);
    }
    template <typename Value>
    void gather(Value*       values,
                Value*       scratch,
                const uint*  indices,
                uint         count,
                uint         blocks,
                uint         threads,
                cudaStream_t stream)
    {
        gatherKernel<<<blocks, threads, 0, stream>>>(
            IndexedGatherView<Value>{values, scratch, indices, count});
        cudaErrCheck(cudaGetLastError());
    }
    template <typename Value>
    void copyBack(Value* values, const Value* scratch, uint count, cudaStream_t stream)
    {
        cudaErrCheck(cudaMemcpyAsync(values,
                                     scratch,
                                     count * sizeof(Value),
                                     cudaMemcpyDeviceToDevice,
                                     stream));
    }
}
namespace grains::cuda
{
    template <typename T>
    struct ParticleSortState
    {
        GrainsMemBuffer<Cells<T, CellOrdering::MORTON>*, MemType::DEVICE> cells;
        grains::detail::ParticleSortBuffers<T, MemType::DEVICE>           buffers;
        GrainsMemBuffer<std::uint64_t, MemType::DEVICE>                   sortedCodes;
        GrainsMemBuffer<uint, MemType::DEVICE>                            sortedIndices;
        GrainsMemBuffer<std::uint8_t, MemType::DEVICE>                    scratch;
        std::array<cudaStream_t, 7>                                       streams{};

        ~ParticleSortState()
        {
            for(auto stream : streams)
                if(stream != nullptr)
                    cudaErrCheck(cudaStreamDestroy(stream));
            // CellsFactory uses device new; cudaFree cannot release that object.
            if(cells.getSize() > 0)
            {
                deleteCells<T><<<1, 1>>>(cells.getData());
                cudaErrCheck(cudaGetLastError());
                cudaErrCheck(cudaDeviceSynchronize());
            }
        }
        void grow(uint count)
        {
            if(count <= buffers.codes.getSize())
                return;
            buffers.grow(count);
            sortedCodes.resize(count);
            sortedIndices.resize(count);
            const auto bytes = detail::particleSortStorageBytes(count);
            if(bytes > scratch.getSize())
                scratch.resize(bytes);
        }
    };
    void ParticleSortStateDeleter::operator()(ParticleSortState<float>* state) const
    {
        delete state;
    }
    void ParticleSortStateDeleter::operator()(ParticleSortState<double>* state) const
    {
        delete state;
    }
    template <typename T>
    ParticleSortStatePtr<T> createParticleSortState(const Vector3<T>& minCorner,
                                                    const Vector3<T>& maxCorner,
                                                    T                 cellSize,
                                                    uint              capacity)
    {
        ParticleSortStatePtr<T> state(new ParticleSortState<T>);
        CellsFactory<T, CellOrdering::MORTON>::template create<MemType::DEVICE>(minCorner,
                                                                                maxCorner,
                                                                                cellSize,
                                                                                state->cells);
        state->grow(capacity);
        for(auto& stream : state->streams)
            cudaErrCheck(cudaStreamCreate(&stream));
        return state;
    }
    template <typename T>
    void sortParticles(ParticleSortState<T>& state, ParticleSortView<T> view)
    {
        if(view.count == 0)
            return;
        state.grow(view.count);
        auto& b = state.buffers;
        uint  blocks, threads;
        computeOptimalThreadsAndBlocks(view.count, GrainsParameters<T>::m_GPU, blocks, threads);
        mortonKernel<<<blocks, threads>>>(ParticleMortonView<T>{view.position,
                                                                b.codes.getData(),
                                                                b.indices.getData(),
                                                                view.count},
                                          state.cells.getData());
        cudaErrCheck(cudaGetLastError());
        cudaErrCheck(cudaDeviceSynchronize());
        detail::sortParticleIndices(b.codes.getData(),
                                    state.sortedCodes.getData(),
                                    b.indices.getData(),
                                    state.sortedIndices.getData(),
                                    view.count,
                                    state.scratch.getData(),
                                    state.scratch.getSize());
        const auto* indices = state.sortedIndices.getData();
        const auto& streams = state.streams;
        gather(view.position,
               b.position.getData(),
               indices,
               view.count,
               blocks,
               threads,
               streams[0]);
        gather(view.velocity,
               b.velocity.getData(),
               indices,
               view.count,
               blocks,
               threads,
               streams[1]);
        gather(view.quaternion,
               b.quaternion.getData(),
               indices,
               view.count,
               blocks,
               threads,
               streams[2]);
        gather(view.torce, b.torce.getData(), indices, view.count, blocks, threads, streams[3]);
        gather(view.bodyTag, b.bodyTag.getData(), indices, view.count, blocks, threads, streams[4]);
        gather(view.localPos,
               b.localPos.getData(),
               indices,
               view.count,
               blocks,
               threads,
               streams[5]);
        gather(view.localQuat,
               b.localQuat.getData(),
               indices,
               view.count,
               blocks,
               threads,
               streams[6]);
        for(auto stream : streams)
            cudaErrCheck(cudaStreamSynchronize(stream));
        copyBack(view.position, b.position.getData(), view.count, streams[0]);
        copyBack(view.velocity, b.velocity.getData(), view.count, streams[1]);
        copyBack(view.quaternion, b.quaternion.getData(), view.count, streams[2]);
        copyBack(view.torce, b.torce.getData(), view.count, streams[3]);
        copyBack(view.bodyTag, b.bodyTag.getData(), view.count, streams[4]);
        copyBack(view.localPos, b.localPos.getData(), view.count, streams[5]);
        copyBack(view.localQuat, b.localQuat.getData(), view.count, streams[6]);
        for(auto stream : streams)
            cudaErrCheck(cudaStreamSynchronize(stream));
    }
    template <typename T>
    void rebuildMasterSlots(MasterSlotView view)
    {
        if(view.count == 0)
            return;
        uint blocks, threads;
        computeOptimalThreadsAndBlocks(view.count, GrainsParameters<T>::m_GPU, blocks, threads);
        masterSlotKernel<<<blocks, threads>>>(view);
        cudaErrCheck(cudaGetLastError());
    }

    template ParticleSortStatePtr<float>
        createParticleSortState<float>(const Vector3<float>&, const Vector3<float>&, float, uint);
    template ParticleSortStatePtr<double> createParticleSortState<double>(const Vector3<double>&,
                                                                          const Vector3<double>&,
                                                                          double,
                                                                          uint);
    template void sortParticles<float>(ParticleSortState<float>&, ParticleSortView<float>);
    template void sortParticles<double>(ParticleSortState<double>&, ParticleSortView<double>);
    template void rebuildMasterSlots<float>(MasterSlotView);
    template void rebuildMasterSlots<double>(MasterSlotView);
}
