#ifndef GRAINS_COLLISION_FILTERING_HH
#define GRAINS_COLLISION_FILTERING_HH

#include "BodyTag.hh"
#include "CollisionDetectionCommon.hh"

struct BVFilterOptions
{
    BoundingVolumeType boundingVolume;
    bool               prebuiltShapes;
};

__HOSTDEVICE__ inline bool
    isIntraCompositePair(IndexPair pair, const uint* tags, uint numComposites)
{
    if(numComposites == 0)
        return false;
    const uint a = tags[pair.x];
    const uint b = tags[pair.y];
    return isSubBody(a) && isSubBody(b) && getCompositeIdx(a) == getCompositeIdx(b);
}

// Non-owning candidate arrays. Filtering examines original pair slots, before
// compaction. Only rejected pairs receive a no-contact sentinel.
template <typename T>
struct BVPairView
{
    const IndexPair*           pairs;
    const RigidBody<T>* const* rigidBodies;
    const BVData<T>*           bounds;
    const uint*                bodyTags;
    ContactInfo<T>*            contacts;
    uint                       count;
    uint                       numComposites;

    __HOSTDEVICE__ bool acceptsComposite(uint i) const
    {
        return !isIntraCompositePair(pairs[i], bodyTags, numComposites);
    }
    __HOSTDEVICE__ bool record(uint i, bool accepted) const
    {
        if(!accepted)
            contacts[i].setOverlapDistance(T(1));
        return accepted;
    }
    template <BoundingVolumeType Bounds, bool Prebuilt>
    __HOSTDEVICE__ bool passes(uint                 i,
                               const Vector3<T>&    relativePosition,
                               const Quaternion<T>& relativeOrientation) const
    {
        const IndexPair pair = pairs[i];
        if constexpr(Prebuilt)
            return filterPairBV_common<T, Bounds>(bounds[getShapeId(bodyTags[pair.x])],
                                                  bounds[getShapeId(bodyTags[pair.y])],
                                                  relativePosition,
                                                  relativeOrientation);
        else
            return filterPairBV_common<T, Bounds>(*rigidBodies[pair.x],
                                                  *rigidBodies[pair.y],
                                                  relativePosition,
                                                  relativeOrientation);
    }
};

template <typename T>
struct RelativeBVView
{
    BVPairView<T>        data;
    const Vector3<T>*    relativePositions;
    const Quaternion<T>* relativeOrientations;

    template <BoundingVolumeType Bounds, bool Prebuilt>
    __HOSTDEVICE__ bool filter(uint i) const
    {
        if(i >= data.count)
            return false;
        if(!data.acceptsComposite(i))
            return data.record(i, false);
        if constexpr(Bounds == BoundingVolumeType::OFF)
            return true;  // Composite-only mode must not read geometry or transforms.
        else
            return data.record(i,
                               data.template passes<Bounds, Prebuilt>(i,
                                                                      relativePositions[i],
                                                                      relativeOrientations[i]));
    }
};

template <typename T>
struct GlobalBVView
{
    BVPairView<T>        data;
    const Vector3<T>*    positions;
    const Quaternion<T>* orientations;

    template <BoundingVolumeType Bounds, bool Prebuilt>
    __HOSTDEVICE__ bool filter(uint i) const
    {
        if(i >= data.count)
            return false;
        if(!data.acceptsComposite(i))
            return data.record(i, false);
        if constexpr(Bounds == BoundingVolumeType::OFF)
            return true;
        else
        {
            const IndexPair  pair             = data.pairs[i];
            const Vector3<T> relativePosition = orientations[pair.x]
                                                << (positions[pair.y] - positions[pair.x]);
            const Quaternion<T> relativeOrientation
                = inverse(orientations[pair.x]) * orientations[pair.y];
            return data.record(
                i,
                data.template passes<Bounds, Prebuilt>(i, relativePosition, relativeOrientation));
        }
    }
};
#endif
