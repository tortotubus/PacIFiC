#ifndef GRAINS_FORCE_COMPUTATION_CUDA_HH
#define GRAINS_FORCE_COMPUTATION_CUDA_HH

#include <cstddef>
#include <memory>

template <typename T>
struct ContactForceView;
template <typename T>
struct ExternalForceView;
template <typename T>
struct CompositeForceView;
class GrainsForceTimer;

namespace grains::cuda
{
    // Defined only in ForceComputationCUDA.cu. Device allocations never appear in
    // the module header, and destruction runs in the CUDA translation unit.
    template <typename T>
    struct ForceState;
    struct ForceStateDeleter
    {
        void operator()(ForceState<float>* state) const;
        void operator()(ForceState<double>* state) const;
    };
    template <typename T>
    using ForceStatePtr = std::unique_ptr<ForceState<T>, ForceStateDeleter>;

    template <typename T>
    ForceStatePtr<T> createForceState(std::size_t capacity);
    template <typename T>
    void resizeForceState(ForceState<T>& state, std::size_t capacity);
    // Device arrays. Work uses the default stream; compaction synchronizes before
    // reading its selected-count result. Other synchronization follows the timer.
    template <typename T>
    void computeContactForces(ForceState<T>&      state,
                              ContactForceView<T> view,
                              bool                compact,
                              GrainsForceTimer&   timer);
    void addExternalForces(ExternalForceView<float> view);
    void addExternalForces(ExternalForceView<double> view);
    void assembleCompositeForces(CompositeForceView<float> view);
    void assembleCompositeForces(CompositeForceView<double> view);
}
#endif
