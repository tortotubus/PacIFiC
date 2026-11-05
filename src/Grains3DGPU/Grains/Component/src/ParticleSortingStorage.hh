#ifndef GRAINS_PARTICLE_SORTING_STORAGE_HH
#define GRAINS_PARTICLE_SORTING_STORAGE_HH

#include "GrainsMemBuffer.hh"
#include "ParticleSorting.hh"

namespace grains::detail
{
    // Private storage shared by the CPU/CUDA state implementations.
    template <typename T, MemType M>
    struct ParticleSortBuffers
    {
        GrainsMemBuffer<std::uint64_t, M> codes;
        GrainsMemBuffer<uint, M>          indices;
        GrainsMemBuffer<Vector3<T>, M>    position;
        GrainsMemBuffer<Kinematics<T>, M> velocity;
        GrainsMemBuffer<Quaternion<T>, M> quaternion;
        GrainsMemBuffer<Torce<T>, M>      torce;
        GrainsMemBuffer<uint, M>          bodyTag;
        GrainsMemBuffer<Vector3<T>, M>    localPos;
        GrainsMemBuffer<Quaternion<T>, M> localQuat;

        void grow(uint count)
        {
            if(count <= codes.getSize())
                return;
            codes.resize(count);
            indices.resize(count);
            position.resize(count);
            velocity.resize(count);
            quaternion.resize(count);
            torce.resize(count);
            bodyTag.resize(count);
            localPos.resize(count);
            localQuat.resize(count);
        }
    };
}
#endif
