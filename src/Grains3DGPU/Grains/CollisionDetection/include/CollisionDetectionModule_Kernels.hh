#ifndef _COLLISIONDETECTIONMODULE_KERNELS_HH_
#define _COLLISIONDETECTIONMODULE_KERNELS_HH_

#include <cstdint>
#include <cuda_runtime.h>

#include "BodyTag.hh"
#include "ContactInfo.hh"
#include "GJK.hh"
#include "GJK_ShapeData.hh"
#include "GrainsParameters.hh"
#include "Quaternion.hh"
#include "RigidBody.hh"
#include "Transform3.hh"
#include "Vector3.hh"

// =================================================================================================
/** @brief CUDA kernel declarations for CollisionDetectionModule.

    Definitions live in CollisionDetectionModule_Kernels.cpp (compiled by nvcc).

    @author A.Yazdani - 2026 - Construction */
// =================================================================================================
/** @name CollisionDetectionModule_Kernels */
//@{

/** @brief BV pre-filter kernel using world-frame positions/quaternions directly.
    Computes pair-relative transforms on-the-fly from world-frame inputs.
    @param rigidBodies   Rigid body array
    @param pairList      List of pairs (full NL)
    @param bodyTags      Per-component body tags (composite membership test)
    @param numComposites Number of composite bodies (0 = no composites, skip the check)
    @param positions     World-frame particle positions
    @param quaternions   World-frame particle quaternions
    @param bvPassFlags   Output pass/fail flags (0 = reject, 1 = pass)
    @param nPairs        Number of pairs */
template <typename T, BoundingVolumeType BVType = BoundingVolumeType::OBB>
__GLOBAL__ void filterPairsBV_Kernel(const RigidBody<T>* const* __RESTRICT__ rigidBodies,
                                     const uint2* __RESTRICT__               pairList,
                                     const uint* __RESTRICT__                bodyTags,
                                     uint                                    numComposites,
                                     const Vector3<T>* __RESTRICT__          positions,
                                     const Quaternion<T>* __RESTRICT__       quaternions,
                                     uint8_t* __RESTRICT__                   bvPassFlags,
                                     const uint                              nPairs);

/** @brief BV pre-filter kernel using pre-built BVData and world-frame positions/quaternions.
    Circumscribed radii are stored in BVData; no RigidBody pointer needed.
    @param bvData        Pre-built BVData array (shapeId-indexed, size = nUniqueShapes)
    @param pairList      List of pairs (full NL)
    @param bodyTags      Per-component body tags (shapeId + composite membership)
    @param numComposites Number of composite bodies
    @param positions     World-frame positions
    @param quaternions   World-frame quaternions
    @param bvPassFlags   Output pass/fail flags
    @param nPairs        Number of pairs */
template <typename T, BoundingVolumeType BVType = BoundingVolumeType::OBB>
__GLOBAL__ void filterPairsBV_Kernel(const BVData<T>* __RESTRICT__     bvData,
                                     const uint2* __RESTRICT__         pairList,
                                     const uint* __RESTRICT__          bodyTags,
                                     uint                              numComposites,
                                     const Vector3<T>* __RESTRICT__    positions,
                                     const Quaternion<T>* __RESTRICT__ quaternions,
                                     uint8_t* __RESTRICT__             bvPassFlags,
                                     const uint                        nPairs);

/** @brief Narrow-phase GJK detection using absolute world-frame positions and quaternions.
    When @p activePairIndices is non-null each thread resolves its original pair index through the
    indirection table and writes compactly; when null the thread ID maps directly (BV-OFF path).
    @param rigidBody         Rigid body array
    @param pairList          Full NL pair list
    @param activePairIndices Compacted index table from CUB, or nullptr to process all pairs
    @param position          World-frame positions
    @param quaternion        World-frame quaternions
    @param contactInfo       Output contact information; compact when activePairIndices != null
    @param compactPairListOut Compact pair output; null on BV-OFF sequential path
    @param nPairs            Number of pairs to process (active or total) */
template <typename T,
          GJKType            GJKVARIANT = GJKType::JOHNSON,
          bool               GJKACC     = false,
          BoundingVolumeType BVType     = BoundingVolumeType::OFF>
__GLOBAL__ void detectCollisionsComponents_Kernel(const RigidBody<T>* const* __RESTRICT__ rigidBody,
                                                  const uint2* __RESTRICT__               pairList,
                                                  const uint* __RESTRICT__       activePairIndices,
                                                  const Vector3<T>* __RESTRICT__ position,
                                                  const Quaternion<T>* __RESTRICT__ quaternion,
                                                  ContactInfo<T>* __RESTRICT__      contactInfo,
                                                  uint2* __RESTRICT__ compactPairListOut,
                                                  const uint          nPairs);
/** @brief Rebuilds the per-composite masterSlot lookup after a Morton sort.
    Only threads whose body tag identifies a sub-body with local index 0 write to the table.
    @param masterSlot  Per-composite master slot array (size = numComposites)
    @param bodyTag     Per-component body tag array (size = nComponents)
    @param nComponents Total number of components (obstacles + particles) */
__GLOBAL__ void rebuildMasterSlot_Kernel(uint* __RESTRICT__       masterSlot,
                                         const uint* __RESTRICT__ bodyTag,
                                         const uint               nComponents);

/** @brief Fills the compact ShapeData table from the rigid body pointer array on the GPU.
    One thread per unique shape; repSlots[k] gives the representative component slot for shape k.
    @param shapeData      Output ShapeData array (size = nUniqueShapes)
    @param rigidBody      Rigid body pointer array (slot-indexed, size >= nComponents)
    @param repSlots       Representative slot for each unique shape (size = nUniqueShapes)
    @param nUniqueShapes  Number of unique shapes (table size) */
template <typename T>
__GLOBAL__ void fillShapeData_Kernel(ShapeData<T>* __RESTRICT__              shapeData,
                                     const RigidBody<T>* const* __RESTRICT__ rigidBody,
                                     const uint* __RESTRICT__                repSlots,
                                     const uint                              nUniqueShapes);

/** @brief Fills the compact BVData table from the rigid body pointer array on the GPU.
    One thread per unique shape; repSlots[k] gives the representative component slot for shape k.
    @param bvData         Output BVData array (size = nUniqueShapes)
    @param rigidBody      Rigid body pointer array (slot-indexed, size >= nComponents)
    @param repSlots       Representative slot for each unique shape (size = nUniqueShapes)
    @param nUniqueShapes  Number of unique shapes (table size) */
template <typename T>
__GLOBAL__ void fillBVData_Kernel(BVData<T>* __RESTRICT__                 bvData,
                                  const RigidBody<T>* const* __RESTRICT__ rigidBody,
                                  const uint* __RESTRICT__                repSlots,
                                  const uint                              nUniqueShapes);

/** @brief Narrow-phase GJK detection with vtable-free support evaluation via ShapeData,
    using absolute world-frame positions and quaternions.
    ShapeData indexed by shapeId (via bodyTags) for deduplication -- no RigidBody needed.
    When activePairIndices is non-null each thread resolves its original pair index via the
    indirection table and writes compactly; when null the thread ID maps directly (BV-OFF).
    @param shapeData         Pre-built ShapeData array (shapeId-indexed, size = nUniqueShapes)
    @param pairList          Full NL pair list
    @param bodyTags          Per-component body tags (encodes shapeId)
    @param activePairIndices Compacted index table from CUB, or nullptr to process all pairs
    @param position          World-frame positions
    @param quaternion        World-frame quaternions
    @param contactInfo       Output contact information; compact when activePairIndices != null
    @param compactPairListOut Compact pair output; null on BV-OFF sequential path
    @param nPairs            Number of pairs to process (active or total) */
template <typename T, GJKType GJKVARIANT = GJKType::JOHNSON, bool GJKACC = false>
__GLOBAL__ __launch_bounds__(256, 2) void detectCollisionsComponents_Kernel(
    const ShapeData<T>* __RESTRICT__  shapeData,
    const uint2* __RESTRICT__         pairList,
    const uint* __RESTRICT__          bodyTags,
    const uint* __RESTRICT__          activePairIndices,
    const Vector3<T>* __RESTRICT__    position,
    const Quaternion<T>* __RESTRICT__ quaternion,
    ContactInfo<T>* __RESTRICT__      contactInfo,
    uint2* __RESTRICT__               compactPairListOut,
    const uint                        nPairs);
//@}

#endif
