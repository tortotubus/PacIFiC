#include "ForceComputationCUDA.hh"
#include "ForceComputation.hh"
#include "ForceCompactionCUDA.hh"
#include "ContactTableDeviceOps.cuh"
#include "GrainsCUDAUtils.hh"
#include "GrainsParameters.hh"
#include "TorceDeviceOps.cuh"

namespace
{
    template <typename T>
    __global__ void computeKernel(PairForceView<T> view, const uint* active, uint count)
    {
        const uint i = blockIdx.x * blockDim.x + threadIdx.x;
        if(i < count)
            view.template compute<grains::cuda::ContactTableOps>(active ? active[i] : i);
    }
    template <typename T>
    __global__ void reduceKernel(PairForceView<T> view, const uint* active, uint count)
    {
        const uint i = blockIdx.x * blockDim.x + threadIdx.x;
        if(i >= count)
            return;
        const uint      pairID = active ? active[i] : i;
        const IndexPair pair   = view.contact.pairs[pairID];
        grains::cuda::TorceOps::addTorceAtomic(view.contact.torces[pair.x], view.first[pairID]);
        grains::cuda::TorceOps::addTorceAtomic(view.contact.torces[pair.y], view.second[pairID]);
    }
    template <typename T>
    __global__ void externalKernel(ExternalForceView<T> view)
    {
        const uint i = blockIdx.x * blockDim.x + threadIdx.x;
        if(i < view.end - view.begin)
            view.add(view.begin + i);
    }
    template <typename T>
    __global__ void compositeKernel(CompositeForceView<T> view)
    {
        const uint offset = blockIdx.x * blockDim.x + threadIdx.x;
        if(offset >= view.end - view.begin)
            return;
        const uint i = view.begin + offset;
        uint       master;
        Vector3<T> force, torque;
        if(view.contribution(i, master, force, torque))
        {
            grains::cuda::TorceOps::addForceAtomic(view.torces[master], force);
            grains::cuda::TorceOps::addTorqueAtomic(view.torces[master], torque);
            view.torces[i].reset();
        }
    }
    template <typename T>
    void add(ExternalForceView<T> view)
    {
        if(view.begin >= view.end)
            return;
        uint blocks, threads;
        computeOptimalThreadsAndBlocks(view.end - view.begin,
                                       GrainsParameters<T>::m_GPU,
                                       blocks,
                                       threads);
        externalKernel<<<blocks, threads>>>(view);
        cudaErrCheck(cudaGetLastError());
    }
    template <typename T>
    void assemble(CompositeForceView<T> view)
    {
        if(view.begin >= view.end)
            return;
        uint blocks, threads;
        computeOptimalThreadsAndBlocks(view.end - view.begin,
                                       GrainsParameters<T>::m_GPU,
                                       blocks,
                                       threads);
        compositeKernel<<<blocks, threads>>>(view);
        cudaErrCheck(cudaGetLastError());
    }
}

namespace grains::cuda
{
    template <typename T>
    struct ForceState
    {
        GrainsMemBuffer<uint, MemType::DEVICE>     active;
        GrainsMemBuffer<uint8_t, MemType::DEVICE>  scratch;
        GrainsMemBuffer<uint, MemType::MAPPED>     selectedCount;
        GrainsMemBuffer<Torce<T>, MemType::DEVICE> first;
        GrainsMemBuffer<Torce<T>, MemType::DEVICE> second;
    };

    void ForceStateDeleter::operator()(ForceState<float>* state) const
    {
        delete state;
    }
    void ForceStateDeleter::operator()(ForceState<double>* state) const
    {
        delete state;
    }

    template <typename T>
    void resizeForceState(ForceState<T>& state, std::size_t capacity)
    {
        state.active.resize(capacity);
        state.scratch.resize(detail::compactionStorageBytes<T>(static_cast<uint>(capacity)));
    }
    template <typename T>
    ForceStatePtr<T> createForceState(std::size_t capacity)
    {
        ForceStatePtr<T> state(new ForceState<T>);
        state->selectedCount.initialize(1);
        *state->selectedCount.getData() = 0;
        resizeForceState(*state, capacity);
        return state;
    }
    template <typename T>
    void computeContactForces(ForceState<T>&      state,
                              ContactForceView<T> view,
                              bool                compact,
                              GrainsForceTimer&   timer)
    {
        uint        count  = view.count;
        const uint* active = nullptr;
        if(compact)
        {
            timer.start(FMStage::FlagAndCompact);
            count = detail::compactContacts(view.contacts,
                                            view.count,
                                            state.active.getData(),
                                            state.selectedCount.getDeviceData(),
                                            state.selectedCount.getData(),
                                            state.scratch.getData(),
                                            state.scratch.getSize());
            timer.stop(FMStage::FlagAndCompact);
            active = state.active.getData();
        }
        if(count == 0)
            return;
        state.first.resize(view.count);
        state.second.resize(view.count);
        const PairForceView<T> pairs{view, state.first.getData(), state.second.getData()};
        uint                   blocks, threads;
        computeOptimalThreadsAndBlocks(count, GrainsParameters<T>::m_GPU, blocks, threads);
        timer.start(FMStage::ComputeForces);
        computeKernel<<<blocks, threads>>>(pairs, active, count);
        cudaErrCheck(cudaGetLastError());
        timer.stop(FMStage::ComputeForces);
        timer.start(FMStage::ReduceTorces);
        reduceKernel<<<blocks, threads>>>(pairs, active, count);
        cudaErrCheck(cudaGetLastError());
        timer.stop(FMStage::ReduceTorces);
    }
    void addExternalForces(ExternalForceView<float> view)
    {
        add(view);
    }
    void addExternalForces(ExternalForceView<double> view)
    {
        add(view);
    }
    void assembleCompositeForces(CompositeForceView<float> view)
    {
        assemble(view);
    }
    void assembleCompositeForces(CompositeForceView<double> view)
    {
        assemble(view);
    }

    template ForceStatePtr<float>  createForceState<float>(std::size_t);
    template ForceStatePtr<double> createForceState<double>(std::size_t);
    template void                  resizeForceState<float>(ForceState<float>&, std::size_t);
    template void                  resizeForceState<double>(ForceState<double>&, std::size_t);
    template void                  computeContactForces<float>(ForceState<float>&,
                                              ContactForceView<float>,
                                              bool,
                                              GrainsForceTimer&);
    template void                  computeContactForces<double>(ForceState<double>&,
                                               ContactForceView<double>,
                                               bool,
                                               GrainsForceTimer&);
}
