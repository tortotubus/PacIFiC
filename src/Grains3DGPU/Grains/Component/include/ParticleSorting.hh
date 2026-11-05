#ifndef GRAINS_PARTICLE_SORTING_HH
#define GRAINS_PARTICLE_SORTING_HH

#include "BodyTag.hh"
#include "Kinematics.hh"
#include "Quaternion.hh"
#include "Torce.hh"
#include <cstdint>

// Arrays point to the particle range only; callers skip the obstacle prefix.
template <typename T>
struct ParticleSortView
{
    Vector3<T>*    position;
    Kinematics<T>* velocity;
    Quaternion<T>* quaternion;
    Torce<T>*      torce;
    uint*          bodyTag;
    Vector3<T>*    localPos;
    Quaternion<T>* localQuat;
    uint           count;
};

template <typename T>
struct ParticleMortonView
{
    const Vector3<T>* position;
    std::uint64_t*    codes;
    uint*             indices;
    uint              count;

    template <typename CellsType>
    __HOSTDEVICE__ void compute(uint i, const CellsType& cells) const
    {
        if(i < count)
        {
            codes[i]   = cells.computeCellHash(position[i]);
            indices[i] = i;
        }
    }
};

// Source and destination must not alias. The same gather is used for each of
// the seven arrays on CPU and on the CUDA backend's independent streams.
template <typename Value>
struct IndexedGatherView
{
    const Value* source;
    Value*       destination;
    const uint*  indices;
    uint         count;

    __HOSTDEVICE__ void gather(uint i) const
    {
        if(i < count)
            destination[i] = source[indices[i]];
    }
};

struct MasterSlotView
{
    const uint* bodyTags;
    uint*       masterSlots;
    uint        count;

    __HOSTDEVICE__ void rebuild(uint i) const
    {
        if(i >= count)
            return;
        const uint tag = bodyTags[i];
        if(isSubBody(tag) && getSubBodyLocalIdx(tag) == 0u)
            masterSlots[getCompositeIdx(tag)] = i;
    }
};
#endif
