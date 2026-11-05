#ifndef GRAINS_PARTICLE_SORTING_CUDA_HH
#define GRAINS_PARTICLE_SORTING_CUDA_HH

#include <memory>

template <typename T>
class Vector3;
template <typename T>
struct ParticleSortView;
struct MasterSlotView;

namespace grains::cuda
{
    template <typename T>
    struct ParticleSortState;
    struct ParticleSortStateDeleter
    {
        void operator()(ParticleSortState<float>* state) const;
        void operator()(ParticleSortState<double>* state) const;
    };
    template <typename T>
    using ParticleSortStatePtr = std::unique_ptr<ParticleSortState<T>, ParticleSortStateDeleter>;

    template <typename T>
    ParticleSortStatePtr<T> createParticleSortState(const Vector3<T>& minCorner,
                                                    const Vector3<T>& maxCorner,
                                                    T                 cellSize,
                                                    unsigned int      capacity);
    // View arrays are device pointers. State owns cells, sort/gather buffers,
    // CUB scratch and seven streams. Growth is automatic; empty work is a no-op.
    // Sorting and all gather/copy streams complete before returning. CUB retains
    // its stable radix ordering for equal Morton keys.
    template <typename T>
    void sortParticles(ParticleSortState<T>& state, ParticleSortView<T> view);
    // Enqueues on the default stream, matching the module's previous behavior.
    template <typename T>
    void rebuildMasterSlots(MasterSlotView view);
}
#endif
