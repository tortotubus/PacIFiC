#ifndef GRAINS_COLLISION_SHAPE_DATA_HH
#define GRAINS_COLLISION_SHAPE_DATA_HH

#include "GJK_ShapeData.hh"

// Borrowed output tables and bodies in one memory space. Representatives map
// shape IDs to component slots; repeated shapes share their first slot.
template <typename T>
struct ShapeTableView
{
    ShapeData<T>*              shapes;
    BVData<T>*                 bounds;
    const RigidBody<T>* const* bodies;
    uint                       count;

    __HOSTDEVICE__ void fillShape(uint i, const uint* representatives) const
    {
        if(i < count)
            fillShapeData(shapes[i], bodies[representatives[i]]);
    }
    __HOSTDEVICE__ void fillBounds(uint i, const uint* representatives) const
    {
        if(i < count)
            fillBVData(bounds[i], bodies[representatives[i]]);
    }
};
#endif
