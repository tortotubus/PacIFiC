#include "CollisionDetectionModule_Kernels.hh"
#include "CollisionDetectionCommon.hh"
#include <cstdint>

// -------------------------------------------------------------------------------------------------
// BV pre-filter using world-frame positions/quaternions; computes relative transforms per-pair
// on-the-fly.
template <typename T, BoundingVolumeType BVType>
__GLOBAL__ void filterPairsBV_Kernel(const RigidBody<T>* const* __RESTRICT__ rigidBodies,
                                     const uint2* __RESTRICT__               pairList,
                                     const uint* __RESTRICT__                bodyTags,
                                     uint                                    numComposites,
                                     const Vector3<T>* __RESTRICT__          positions,
                                     const Quaternion<T>* __RESTRICT__       quaternions,
                                     uint8_t* __RESTRICT__                   bvPassFlags,
                                     const uint                              nPairs)
{
    uint tID = blockIdx.x * blockDim.x + threadIdx.x;
    if(tID >= nPairs)
        return;

    const uint2 pair = pairList[tID];

    // Reject intra-composite pairs (same composite, different sub-bodies)
    if(numComposites > 0)
    {
        const uint tagA = bodyTags[pair.x];
        const uint tagB = bodyTags[pair.y];
        if(isSubBody(tagA) && isSubBody(tagB) && getCompositeIdx(tagA) == getCompositeIdx(tagB))
        {
            bvPassFlags[tID] = 0;
            return;
        }
    }

    if constexpr(BVType == BoundingVolumeType::OFF)
    {
        bvPassFlags[tID] = 1;
    }
    else
    {
        const RigidBody<T>& rbA   = *(rigidBodies[pair.x]);
        const RigidBody<T>& rbB   = *(rigidBodies[pair.y]);
        const Vector3<T>    v_b2a = quaternions[pair.x] << (positions[pair.y] - positions[pair.x]);
        const Quaternion<T> q_b2a = inverse(quaternions[pair.x]) * quaternions[pair.y];
        bvPassFlags[tID]          = filterPairBV_common<T, BVType>(rbA, rbB, v_b2a, q_b2a) ? 1 : 0;
    }
}

// -------------------------------------------------------------------------------------------------
// Narrow-phase GJK detection using absolute world-frame positions and quaternions.
// When activePairIndices is non-null, writes compactly (contactInfo[tID], compactPairListOut[tID]);
// when null, writes sequentially (BV-off path, all pairs, tID==pairIdx).
template <typename T, GJKType GJKVARIANT, bool GJKACC, BoundingVolumeType BVType>
__GLOBAL__ void detectCollisionsComponents_Kernel(const RigidBody<T>* const* __RESTRICT__ rigidBody,
                                                  const uint2* __RESTRICT__               pairList,
                                                  const uint* __RESTRICT__       activePairIndices,
                                                  const Vector3<T>* __RESTRICT__ position,
                                                  const Quaternion<T>* __RESTRICT__ quaternion,
                                                  ContactInfo<T>* __RESTRICT__      contactInfo,
                                                  uint2* __RESTRICT__ compactPairListOut,
                                                  const uint          nPairs)
{
    uint tID = blockIdx.x * blockDim.x + threadIdx.x;
    if(tID >= nPairs)
        return;

    const uint pairIdx  = (activePairIndices != nullptr) ? activePairIndices[tID] : tID;
    const uint writeIdx = (compactPairListOut != nullptr) ? tID : pairIdx;
    if(compactPairListOut != nullptr)
        compactPairListOut[tID] = pairList[pairIdx];
    detectCollisionsComponents_common<T, GJKVARIANT, GJKACC, BVType>(pairList,
                                                                     rigidBody,
                                                                     position,
                                                                     quaternion,
                                                                     contactInfo,
                                                                     pairIdx,
                                                                     writeIdx);
}

// -------------------------------------------------------------------------------------------------
// Rebuilds masterSlot lookup after a Morton sort
__GLOBAL__ void rebuildMasterSlot_Kernel(uint* __RESTRICT__       masterSlot,
                                         const uint* __RESTRICT__ bodyTag,
                                         const uint               nComponents)
{
    uint cID = blockIdx.x * blockDim.x + threadIdx.x;
    if(cID >= nComponents)
        return;

    const uint tag = bodyTag[cID];
    if(isSubBody(tag) && getSubBodyLocalIdx(tag) == 0u)
        masterSlot[getCompositeIdx(tag)] = cID;
}

// -------------------------------------------------------------------------------------------------
// BV pre-filter using pre-built BVData and world-frame transforms (vtable-free).
template <typename T, BoundingVolumeType BVType>
__GLOBAL__ void filterPairsBV_Kernel(const BVData<T>* __RESTRICT__     bvData,
                                     const uint2* __RESTRICT__         pairList,
                                     const uint* __RESTRICT__          bodyTags,
                                     uint                              numComposites,
                                     const Vector3<T>* __RESTRICT__    positions,
                                     const Quaternion<T>* __RESTRICT__ quaternions,
                                     uint8_t* __RESTRICT__             bvPassFlags,
                                     const uint                        nPairs)
{
    uint tID = blockIdx.x * blockDim.x + threadIdx.x;
    if(tID >= nPairs)
        return;

    const uint2 pair = pairList[tID];

    if(numComposites > 0)
    {
        const uint tagA = bodyTags[pair.x];
        const uint tagB = bodyTags[pair.y];
        if(isSubBody(tagA) && isSubBody(tagB) && getCompositeIdx(tagA) == getCompositeIdx(tagB))
        {
            bvPassFlags[tID] = 0;
            return;
        }
    }

    if constexpr(BVType == BoundingVolumeType::OFF)
    {
        bvPassFlags[tID] = 1;
    }
    else
    {
        const BVData<T>&    bvA   = bvData[getShapeId(bodyTags[pair.x])];
        const BVData<T>&    bvB   = bvData[getShapeId(bodyTags[pair.y])];
        const Vector3<T>    v_b2a = quaternions[pair.x] << (positions[pair.y] - positions[pair.x]);
        const Quaternion<T> q_b2a = inverse(quaternions[pair.x]) * quaternions[pair.y];
        bvPassFlags[tID]          = filterPairBV_common<T, BVType>(bvA, bvB, v_b2a, q_b2a) ? 1 : 0;
    }
}

// -------------------------------------------------------------------------------------------------
// Fills the compact ShapeData table; one thread per unique shape.
// Thread k reads rigidBody[repSlots[k]] and fills shapeData[k].
template <typename T>
__GLOBAL__ void fillShapeData_Kernel(ShapeData<T>* __RESTRICT__              shapeData,
                                     const RigidBody<T>* const* __RESTRICT__ rigidBody,
                                     const uint* __RESTRICT__                repSlots,
                                     const uint                              nUniqueShapes)
{
    uint tID = blockIdx.x * blockDim.x + threadIdx.x;
    if(tID >= nUniqueShapes)
        return;
    fillShapeData(shapeData[tID], rigidBody[repSlots[tID]]);
}

// -------------------------------------------------------------------------------------------------
// Fills the compact BVData table; one thread per unique shape.
// Thread k reads rigidBody[repSlots[k]] and fills bvData[k].
template <typename T>
__GLOBAL__ void fillBVData_Kernel(BVData<T>* __RESTRICT__                 bvData,
                                  const RigidBody<T>* const* __RESTRICT__ rigidBody,
                                  const uint* __RESTRICT__                repSlots,
                                  const uint                              nUniqueShapes)
{
    uint tID = blockIdx.x * blockDim.x + threadIdx.x;
    if(tID >= nUniqueShapes)
        return;
    fillBVData(bvData[tID], rigidBody[repSlots[tID]]);
}

// -------------------------------------------------------------------------------------------------
// Narrow-phase GJK detection with vtable-free support evaluation via ShapeData (world frame).
// ShapeData indexed by shapeId (via bodyTags); no RigidBody pointer needed.
// When activePairIndices is non-null each thread resolves its pair index through the indirection
// table (BV-compacted path). When null the thread ID is used directly (BV-off path, all pairs).
template <typename T, GJKType GJKVARIANT, bool GJKACC>
__GLOBAL__ __launch_bounds__(256, 2) void detectCollisionsComponents_Kernel(
    const ShapeData<T>* __RESTRICT__  shapeData,
    const uint2* __RESTRICT__         pairList,
    const uint* __RESTRICT__          bodyTags,
    const uint* __RESTRICT__          activePairIndices,
    const Vector3<T>* __RESTRICT__    position,
    const Quaternion<T>* __RESTRICT__ quaternion,
    ContactInfo<T>* __RESTRICT__      contactInfo,
    uint2* __RESTRICT__               compactPairListOut,
    const uint                        nPairs)
{
    uint tID = blockIdx.x * blockDim.x + threadIdx.x;
    if(tID >= nPairs)
        return;

    const uint pairIdx  = (activePairIndices != nullptr) ? activePairIndices[tID] : tID;
    const uint writeIdx = (compactPairListOut != nullptr) ? tID : pairIdx;
    if(compactPairListOut != nullptr)
        compactPairListOut[tID] = pairList[pairIdx];
    detectCollisionsComponents_common<T, GJKVARIANT, GJKACC, BoundingVolumeType::OFF>(pairList,
                                                                                      shapeData,
                                                                                      bodyTags,
                                                                                      position,
                                                                                      quaternion,
                                                                                      contactInfo,
                                                                                      pairIdx,
                                                                                      writeIdx);
}

// -------------------------------------------------------------------------------------------------
// Explicit template instantiations
#define FOR_EACH_FP_TYPE(MACRO) \
    MACRO(float)                \
    MACRO(double)

// 1) Narrow-phase kernels with RigidBody input
#define INSTANTIATE_DETECT_RB(T, GJK, ACC, BV)                                   \
    template __GLOBAL__ void detectCollisionsComponents_Kernel<T, GJK, ACC, BV>( \
        const RigidBody<T>* const*,                                              \
        const uint2*,                                                            \
        const uint*,                                                             \
        const Vector3<T>*,                                                       \
        const Quaternion<T>*,                                                    \
        ContactInfo<T>*,                                                         \
        uint2*,                                                                  \
        const uint);

#define INSTANTIATE_DETECT_RB_FOR_BV(T, BV)                    \
    INSTANTIATE_DETECT_RB(T, GJKType::JOHNSON, false, BV)      \
    INSTANTIATE_DETECT_RB(T, GJKType::JOHNSON, true, BV)       \
    INSTANTIATE_DETECT_RB(T, GJKType::SIGNEDVOLUME, false, BV) \
    INSTANTIATE_DETECT_RB(T, GJKType::SIGNEDVOLUME, true, BV)  \
    INSTANTIATE_DETECT_RB(T, GJKType::SPHERE, false, BV)

#define INSTANTIATE_DETECT_RB_FOR_TYPE(T)                    \
    INSTANTIATE_DETECT_RB_FOR_BV(T, BoundingVolumeType::OFF) \
    INSTANTIATE_DETECT_RB_FOR_BV(T, BoundingVolumeType::OBB) \
    INSTANTIATE_DETECT_RB_FOR_BV(T, BoundingVolumeType::OBC)

FOR_EACH_FP_TYPE(INSTANTIATE_DETECT_RB_FOR_TYPE)

#undef INSTANTIATE_DETECT_RB_FOR_TYPE
#undef INSTANTIATE_DETECT_RB_FOR_BV
#undef INSTANTIATE_DETECT_RB

// 2) BV pre-filter kernels
#define INSTANTIATE_FILTER_RB(T, BV)                                                 \
    template __GLOBAL__ void filterPairsBV_Kernel<T, BV>(const RigidBody<T>* const*, \
                                                         const uint2*,               \
                                                         const uint*,                \
                                                         uint,                       \
                                                         const Vector3<T>*,          \
                                                         const Quaternion<T>*,       \
                                                         uint8_t*,                   \
                                                         const uint);

#define INSTANTIATE_FILTER_RB_FOR_TYPE(T)             \
    INSTANTIATE_FILTER_RB(T, BoundingVolumeType::OBB) \
    INSTANTIATE_FILTER_RB(T, BoundingVolumeType::OBC) \
    INSTANTIATE_FILTER_RB(T, BoundingVolumeType::OFF)

FOR_EACH_FP_TYPE(INSTANTIATE_FILTER_RB_FOR_TYPE)

#undef INSTANTIATE_FILTER_RB_FOR_TYPE
#undef INSTANTIATE_FILTER_RB

#define INSTANTIATE_FILTER_BVDATA(T, BV)                                       \
    template __GLOBAL__ void filterPairsBV_Kernel<T, BV>(const BVData<T>*,     \
                                                         const uint2*,         \
                                                         const uint*,          \
                                                         uint,                 \
                                                         const Vector3<T>*,    \
                                                         const Quaternion<T>*, \
                                                         uint8_t*,             \
                                                         const uint);

#define INSTANTIATE_FILTER_BVDATA_FOR_TYPE(T)             \
    INSTANTIATE_FILTER_BVDATA(T, BoundingVolumeType::OBB) \
    INSTANTIATE_FILTER_BVDATA(T, BoundingVolumeType::OBC) \
    INSTANTIATE_FILTER_BVDATA(T, BoundingVolumeType::OFF)

FOR_EACH_FP_TYPE(INSTANTIATE_FILTER_BVDATA_FOR_TYPE)

#undef INSTANTIATE_FILTER_BVDATA_FOR_TYPE
#undef INSTANTIATE_FILTER_BVDATA

// 3) Shape/BV table build kernels
#define INSTANTIATE_FILL_SHAPE(T)                                                \
    template __GLOBAL__ void fillShapeData_Kernel<T>(ShapeData<T>*,              \
                                                     const RigidBody<T>* const*, \
                                                     const uint*,                \
                                                     const uint);

#define INSTANTIATE_FILL_BV(T)                                                \
    template __GLOBAL__ void fillBVData_Kernel<T>(BVData<T>*,                 \
                                                  const RigidBody<T>* const*, \
                                                  const uint*,                \
                                                  const uint);

FOR_EACH_FP_TYPE(INSTANTIATE_FILL_SHAPE)
FOR_EACH_FP_TYPE(INSTANTIATE_FILL_BV)

#undef INSTANTIATE_FILL_BV
#undef INSTANTIATE_FILL_SHAPE

// 4) Narrow-phase kernels with prebuilt ShapeData input
#define INSTANTIATE_DETECT_PREBUILT(T, GJK, ACC)                                                  \
    template __GLOBAL__ void detectCollisionsComponents_Kernel<T, GJK, ACC>(const ShapeData<T>*,  \
                                                                            const uint2*,         \
                                                                            const uint*,          \
                                                                            const uint*,          \
                                                                            const Vector3<T>*,    \
                                                                            const Quaternion<T>*, \
                                                                            ContactInfo<T>*,      \
                                                                            uint2*,               \
                                                                            const uint);

#define INSTANTIATE_DETECT_PREBUILT_FOR_TYPE(T)                  \
    INSTANTIATE_DETECT_PREBUILT(T, GJKType::JOHNSON, false)      \
    INSTANTIATE_DETECT_PREBUILT(T, GJKType::JOHNSON, true)       \
    INSTANTIATE_DETECT_PREBUILT(T, GJKType::SIGNEDVOLUME, false) \
    INSTANTIATE_DETECT_PREBUILT(T, GJKType::SIGNEDVOLUME, true)  \
    INSTANTIATE_DETECT_PREBUILT(T, GJKType::SPHERE, false)

FOR_EACH_FP_TYPE(INSTANTIATE_DETECT_PREBUILT_FOR_TYPE)

#undef INSTANTIATE_DETECT_PREBUILT_FOR_TYPE
#undef INSTANTIATE_DETECT_PREBUILT
#undef FOR_EACH_FP_TYPE
