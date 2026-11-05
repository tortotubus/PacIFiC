#ifndef GRAINS_PARTICLE_SORTING_CPU_HH
#define GRAINS_PARTICLE_SORTING_CPU_HH

#include <memory>

template <typename T>
class Vector3;
template <typename T>
struct ParticleSortView;
struct MasterSlotView;

namespace grains::cpu
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
    // State owns the cells and reusable scratch arrays. It grows if needed.
    // Equal Morton keys retain the existing CPU std::sort ordering policy.
    template <typename T>
    void sortParticles(ParticleSortState<T>& state, ParticleSortView<T> view);
    void rebuildMasterSlots(MasterSlotView view);
}
#endif
