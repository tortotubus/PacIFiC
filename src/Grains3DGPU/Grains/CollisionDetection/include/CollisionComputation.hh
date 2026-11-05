#ifndef GRAINS_COLLISION_COMPUTATION_HH
#define GRAINS_COLLISION_COMPUTATION_HH

#include "BodyTag.hh"
#include "CollisionDetectionCommon.hh"
#include "CollisionFiltering.hh"

struct NarrowPhaseOptions
{
    NarrowPhaseType    algorithm;
    BoundingVolumeType boundingVolume;
    bool               acceleration;
    bool               prebuiltShapes;
};

template <typename T>
struct RelativeTransformView
{
    const IndexPair*     pairs;
    const Vector3<T>*    positions;
    const Quaternion<T>* orientations;
    Vector3<T>*          relativePositions;
    Quaternion<T>*       relativeOrientations;
    uint                 count;

    __HOSTDEVICE__ void compute(uint i) const
    {
        if(i < count)
            computeRelativeTransformations_common(pairs,
                                                  positions,
                                                  orientations,
                                                  relativePositions,
                                                  relativeOrientations,
                                                  i);
    }
};

// count is the number of work items. activeIndices, when present, maps those
// items to original pair slots; results always go to those original slots.
template <typename T>
struct CollisionPairView
{
    const IndexPair*           pairs;
    const RigidBody<T>* const* rigidBodies;
    const ShapeData<T>*        shapes;
    const uint*                bodyTags;
    ContactInfo<T>*            contacts;
    const uint*                activeIndices;
    uint                       count;
    uint                       numComposites;

    __HOSTDEVICE__ uint pairIndex(uint workIndex) const
    {
        return activeIndices ? activeIndices[workIndex] : workIndex;
    }

    __HOSTDEVICE__ bool accepts(uint pair) const
    {
        if(isIntraCompositePair(pairs[pair], bodyTags, numComposites))
        {
            contacts[pair].setOverlapDistance(T(1));
            return false;
        }
        return true;
    }
};

template <typename T>
struct RelativeCollisionView
{
    CollisionPairView<T> data;
    const Vector3<T>*    relativePositions;
    const Quaternion<T>* relativeOrientations;

    template <GJKType Variant, bool Accelerated, BoundingVolumeType BoundingVolume>
    __HOSTDEVICE__ void detect(uint workIndex, bool prebuilt) const
    {
        if(workIndex >= data.count)
            return;
        const uint pair = data.pairIndex(workIndex);
        if(!data.accepts(pair))
            return;
        if(prebuilt)
            detectCollisionsComponents_common<T, Variant, Accelerated, BoundingVolume>(
                data.pairs,
                data.shapes,
                data.bodyTags,
                relativePositions,
                relativeOrientations,
                data.contacts,
                pair);
        else
            detectCollisionsComponents_common<T, Variant, Accelerated, BoundingVolume>(
                data.pairs,
                data.rigidBodies,
                relativePositions,
                relativeOrientations,
                data.contacts,
                pair);
    }
};

template <typename T>
struct GlobalCollisionView
{
    CollisionPairView<T> data;
    const Vector3<T>*    positions;
    const Quaternion<T>* orientations;

    template <GJKType Variant, bool Accelerated, BoundingVolumeType BoundingVolume>
    __HOSTDEVICE__ void detect(uint workIndex, bool prebuilt) const
    {
        if(workIndex >= data.count)
            return;
        const uint pair = data.pairIndex(workIndex);
        if(!data.accepts(pair))
            return;
        if(prebuilt)
            detectCollisionsComponentsGlobal_common<T, Variant, Accelerated, BoundingVolume>(
                data.pairs,
                data.shapes,
                data.bodyTags,
                positions,
                orientations,
                data.contacts,
                pair);
        else
            detectCollisionsComponentsGlobal_common<T, Variant, Accelerated, BoundingVolume>(
                data.pairs,
                data.rigidBodies,
                positions,
                orientations,
                data.contacts,
                pair);
    }
};

template <typename T>
struct ContactTransformView
{
    const IndexPair*     pairs;
    const Vector3<T>*    positions;
    const Quaternion<T>* orientations;
    ContactInfo<T>*      localContacts;
    ContactInfo<T>*      worldContacts;
    uint                 count;

    __HOSTDEVICE__ void transform(uint i) const
    {
        if(i < count)
            transformContactInfo_common(pairs,
                                        positions,
                                        orientations,
                                        localContacts,
                                        worldContacts,
                                        i);
    }
};
#endif
