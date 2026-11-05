#include "ParticleSortingCPU.hh"
#include "Cells.hh"
#include "ParticleSortingStorage.hh"
#include <algorithm>

namespace grains::cpu
{
    template <typename T>
    struct ParticleSortState
    {
        Cells<T, CellOrdering::MORTON>                        cells;
        grains::detail::ParticleSortBuffers<T, MemType::HOST> buffers;

        ParticleSortState(const Vector3<T>& minCorner, const Vector3<T>& maxCorner, T size)
            : cells(minCorner, maxCorner, size)
        {
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
        GAssert(cellSize > 0, "Cell size must be positive! Aborting Grains!");
        ParticleSortStatePtr<T> state(new ParticleSortState<T>(minCorner, maxCorner, cellSize));
        state->buffers.grow(capacity);
        return state;
    }
    namespace
    {
        template <typename Value>
        void reorder(Value* values, Value* scratch, const uint* indices, uint count)
        {
            const IndexedGatherView<Value> view{values, scratch, indices, count};
            for(uint i = 0; i < count; ++i)
                view.gather(i);
            std::copy_n(scratch, count, values);
        }
    }
    template <typename T>
    void sortParticles(ParticleSortState<T>& state, ParticleSortView<T> view)
    {
        if(view.count == 0)
            return;
        auto& buffers = state.buffers;
        buffers.grow(view.count);
        auto*                       codes   = buffers.codes.getData();
        auto*                       indices = buffers.indices.getData();
        const ParticleMortonView<T> morton{view.position, codes, indices, view.count};
        for(uint i = 0; i < view.count; ++i)
            morton.compute(i, state.cells);
        std::sort(indices, indices + view.count, [codes](uint a, uint b) {
            return codes[a] < codes[b];
        });
        reorder(view.position, buffers.position.getData(), indices, view.count);
        reorder(view.velocity, buffers.velocity.getData(), indices, view.count);
        reorder(view.quaternion, buffers.quaternion.getData(), indices, view.count);
        reorder(view.torce, buffers.torce.getData(), indices, view.count);
        reorder(view.bodyTag, buffers.bodyTag.getData(), indices, view.count);
        reorder(view.localPos, buffers.localPos.getData(), indices, view.count);
        reorder(view.localQuat, buffers.localQuat.getData(), indices, view.count);
    }
    void rebuildMasterSlots(MasterSlotView view)
    {
        for(uint i = 0; i < view.count; ++i)
            view.rebuild(i);
    }

    template ParticleSortStatePtr<float>
        createParticleSortState<float>(const Vector3<float>&, const Vector3<float>&, float, uint);
    template ParticleSortStatePtr<double> createParticleSortState<double>(const Vector3<double>&,
                                                                          const Vector3<double>&,
                                                                          double,
                                                                          uint);
    template void sortParticles<float>(ParticleSortState<float>&, ParticleSortView<float>);
    template void sortParticles<double>(ParticleSortState<double>&, ParticleSortView<double>);
}
