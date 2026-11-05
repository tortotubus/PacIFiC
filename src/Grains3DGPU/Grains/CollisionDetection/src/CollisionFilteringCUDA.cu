#include "CollisionFilteringCUDA.hh"
#include "CollisionFiltering.hh"
#include "CollisionCompactionCUDA.hh"
#include "CollisionDispatch.hh"
#include "GrainsCUDAUtils.hh"
#include "GrainsMemBuffer.hh"
#include "GrainsParameters.hh"

namespace
{
    template <BoundingVolumeType Bounds, bool Prebuilt, typename View>
    __global__ void filterKernel(View view, std::uint8_t* flags)
    {
        const uint i = blockIdx.x * blockDim.x + threadIdx.x;
        if(i < view.data.count)
            flags[i] = static_cast<std::uint8_t>(view.template filter<Bounds, Prebuilt>(i));
    }
}
namespace grains::cuda
{
    template <typename T>
    struct CollisionFilterState
    {
        GrainsMemBuffer<std::uint8_t, MemType::DEVICE> flags;
        GrainsMemBuffer<uint, MemType::DEVICE>         indices;
        GrainsMemBuffer<int, MemType::MAPPED>          selectedCount;
        GrainsMemBuffer<std::uint8_t, MemType::DEVICE> scratch;
    };
    void CollisionFilterStateDeleter::operator()(CollisionFilterState<float>* state) const
    {
        delete state;
    }
    void CollisionFilterStateDeleter::operator()(CollisionFilterState<double>* state) const
    {
        delete state;
    }

    template <typename T>
    void resizeCollisionFilterState(CollisionFilterState<T>& state, std::size_t capacity)
    {
        state.flags.resize(capacity);
        state.indices.resize(capacity);
        const std::size_t bytes
            = detail::collisionCompactionStorageBytes(static_cast<uint>(capacity));
        if(bytes > state.scratch.getSize())
            state.scratch.resize(bytes);
    }
    template <typename T>
    CollisionFilterStatePtr<T> createCollisionFilterState(std::size_t capacity)
    {
        CollisionFilterStatePtr<T> state(new CollisionFilterState<T>);
        state->selectedCount.initialize(1);
        *state->selectedCount.getData() = 0;
        resizeCollisionFilterState(*state, capacity);
        return state;
    }
    namespace
    {
        template <typename T, typename View>
        SelectedPairIndices
            filter(CollisionFilterState<T>& state, View view, BVFilterOptions options)
        {
            if(view.data.count == 0)
            {
                *state.selectedCount.getData() = 0;
                return {state.indices.getData(), 0};
            }
            if(view.data.count > state.flags.getSize())
                resizeCollisionFilterState(state, view.data.count);
            uint blocks, threads;
            computeOptimalThreadsAndBlocks(view.data.count,
                                           GrainsParameters<T>::m_GPU,
                                           blocks,
                                           threads);
            grains::detail::dispatchBoundingVolume(options.boundingVolume, [&](auto bounds) {
                if(options.prebuiltShapes)
                    filterKernel<decltype(bounds)::value, true>
                        <<<blocks, threads>>>(view, state.flags.getData());
                else
                    filterKernel<decltype(bounds)::value, false>
                        <<<blocks, threads>>>(view, state.flags.getData());
            });
            cudaErrCheck(cudaGetLastError());
            // Retain the existing filter-stage synchronization before CUB selection.
            cudaErrCheck(cudaDeviceSynchronize());
            const uint count = detail::compactCollisionFlags(state.flags.getData(),
                                                             view.data.count,
                                                             state.indices.getData(),
                                                             state.selectedCount.getDeviceData(),
                                                             state.selectedCount.getData(),
                                                             state.scratch.getData(),
                                                             state.scratch.getSize());
            return {state.indices.getData(), count};
        }
    }
    template <typename T>
    SelectedPairIndices filterRelativePairs(CollisionFilterState<T>& state,
                                            RelativeBVView<T>        view,
                                            BVFilterOptions          options)
    {
        return filter(state, view, options);
    }
    template <typename T>
    SelectedPairIndices filterGlobalPairs(CollisionFilterState<T>& state,
                                          GlobalBVView<T>          view,
                                          BVFilterOptions          options)
    {
        return filter(state, view, options);
    }

    template CollisionFilterStatePtr<float>  createCollisionFilterState<float>(std::size_t);
    template CollisionFilterStatePtr<double> createCollisionFilterState<double>(std::size_t);
    template void resizeCollisionFilterState<float>(CollisionFilterState<float>&, std::size_t);
    template void resizeCollisionFilterState<double>(CollisionFilterState<double>&, std::size_t);
    template SelectedPairIndices filterRelativePairs<float>(CollisionFilterState<float>&,
                                                            RelativeBVView<float>,
                                                            BVFilterOptions);
    template SelectedPairIndices filterRelativePairs<double>(CollisionFilterState<double>&,
                                                             RelativeBVView<double>,
                                                             BVFilterOptions);
    template SelectedPairIndices filterGlobalPairs<float>(CollisionFilterState<float>&,
                                                          GlobalBVView<float>,
                                                          BVFilterOptions);
    template SelectedPairIndices filterGlobalPairs<double>(CollisionFilterState<double>&,
                                                           GlobalBVView<double>,
                                                           BVFilterOptions);
}
