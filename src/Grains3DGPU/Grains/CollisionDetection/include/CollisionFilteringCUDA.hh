#ifndef GRAINS_COLLISION_FILTERING_CUDA_HH
#define GRAINS_COLLISION_FILTERING_CUDA_HH

#include <cstddef>
#include <memory>

struct BVFilterOptions;
template <typename T>
struct RelativeBVView;
template <typename T>
struct GlobalBVView;

// indices refers to CUDA state-owned device memory, not an owning allocation.
// It is valid until that state is resized, filtered again, or destroyed.
struct SelectedPairIndices
{
    const unsigned int* indices = nullptr;
    unsigned int        count   = 0;
};

namespace grains::cuda
{
    template <typename T>
    struct CollisionFilterState;
    struct CollisionFilterStateDeleter
    {
        void operator()(CollisionFilterState<float>* state) const;
        void operator()(CollisionFilterState<double>* state) const;
    };
    template <typename T>
    using CollisionFilterStatePtr
        = std::unique_ptr<CollisionFilterState<T>, CollisionFilterStateDeleter>;

    template <typename T>
    CollisionFilterStatePtr<T> createCollisionFilterState(std::size_t capacity);
    template <typename T>
    void resizeCollisionFilterState(CollisionFilterState<T>& state, std::size_t capacity);
    // Candidate views refer to device arrays. Filtering and stable compaction run
    // on the default stream and synchronize before returning the selected count.
    // The backend grows its buffers when needed; empty work launches no kernels.
    template <typename T>
    SelectedPairIndices filterRelativePairs(CollisionFilterState<T>& state,
                                            RelativeBVView<T>        view,
                                            BVFilterOptions          options);
    template <typename T>
    SelectedPairIndices filterGlobalPairs(CollisionFilterState<T>& state,
                                          GlobalBVView<T>          view,
                                          BVFilterOptions          options);
}
#endif
