#ifndef _COLLISIONDETECTIONMODULE_HH_
#define _COLLISIONDETECTIONMODULE_HH_

#include <cstdint>
#include <memory>

#include "BodyTag.hh"
#include "CollisionDetectionCommon.hh"
#include "CollisionFilteringCUDA.hh"
#include "ContactInfo.hh"
#include "GJK_ShapeData.hh"
#include "GrainsMemBuffer.hh"
#include "GrainsParameters.hh"
#include "GrainsUtils.hh"
#include "Kinematics.hh"
#include "NeighborList.hh"
#include "NeighborListFactory.hh"
#include "ParticleSorter.hh"
#include "Quaternion.hh"
#include "RigidBody.hh"
#include "Torce.hh"
#include "Vector3.hh"

namespace grains::detail
{
    template <typename T, MemType M>
    struct CollisionFilterStorage
    {
    };
    template <typename T>
    struct CollisionFilterStorage<T, MemType::DEVICE>
    {
        cuda::CollisionFilterStatePtr<T> state;
    };
}

// =================================================================================================
/** @brief The class CollisionDetectionModule.

    Encapsulates the full collision detection pipeline:
      1. Neighbor list update (broad phase)
      2. Relative transformation computation
      3. Optional bounding volume pre-filter (Sphere/OBB/OBC)
      4. Narrow phase GJK (closestPointsRigidBodies)
      5. Contact information world-frame transform

    ComponentManager owns a CollisionDetectionModule by unique_ptr and delegates all detection
    work to it by passing raw non-owning pointers to the particle data arrays.  No particle data
    is copied, and the module reads/writes directly into ComponentManager's buffers.

    @author A.Yazdani - 2026 - Construction */
// =================================================================================================
template <typename T, MemType M = MemType::HOST>
class CollisionDetectionModule
{
private:
    /** @name Parameters */
    //@{
    /** \brief Pre-built ShapeData array for vtable-free GJK support evaluation.
        ShapeId-indexed (one entry per unique shape, size = nUniqueShapes). */
    GrainsMemBuffer<ShapeData<T>, M> m_shapeData;
    /** \brief Pre-built BVData array for vtable-free bounding-volume queries.
        ShapeId-indexed (one entry per unique shape, size = nUniqueShapes). */
    GrainsMemBuffer<BVData<T>, M> m_bvData;
    /** \brief Particle sorter for Morton code-based reordering of particle arrays */
    ParticleSorter<T, M> m_particleSorter;
    /** \brief Neighbor list (broad phase) */
    std::unique_ptr<NeighborList<T, M>> m_neighborList;
    /** \brief Per-pair relative position of B in A-local frame */
    GrainsMemBuffer<Vector3<T>, M> m_relPosition;
    /** \brief Per-pair relative quaternion of B in A-local frame */
    GrainsMemBuffer<Quaternion<T>, M> m_relQuaternion;
    /** \brief Per-pair contact information in A-local frame (output of GJK) */
    GrainsMemBuffer<ContactInfo<T>, M> m_contactInfoLocal;
    /** \brief CUDA-owned filter flags, selected indices, mapped count and scratch.
        HOST modules have no filter state or CUDA resource destruction. */
    grains::detail::CollisionFilterStorage<T, M> m_filterBackend;
    /** \brief Number of unique shapes (= max shapeId+1 from bodyTags). */
    uint m_nUniqueShapes = 0;
    //@}

public:
    /** @name Constructors */
    //@{
    /** @brief Default constructor (forbidden) */
    CollisionDetectionModule() = delete;

    /** @brief Constructor -- builds the NeighborList internally.
        @param rigidBody    Pointer to the rigid body buffer
        @param positions    Position buffer (used for initial cell assignment)
        @param orientations Quaternion buffer
        @param bodyTags     Body-tag buffer (encodes shapeId; used to build compact shape tables)
        @param CD           Collision detection parameters
        @param nObstacles   Number of obstacles
        @param nParticles   Number of moving particles */
    CollisionDetectionModule(const GrainsMemBuffer<RigidBody<T>*, M>* rigidBody,
                             const GrainsMemBuffer<Vector3<T>, M>&    positions,
                             const GrainsMemBuffer<Quaternion<T>, M>& orientations,
                             const GrainsMemBuffer<uint, M>&          bodyTags,
                             const uint*                              hostBodyTags,
                             const CollisionDetectionParameters<T>&   CD,
                             uint                                     nObstacles,
                             uint                                     nParticles);

    /** @brief Destructor */
    ~CollisionDetectionModule() = default;

    /** @brief Deleted copy constructor */
    CollisionDetectionModule(const CollisionDetectionModule&) = delete;

    /** @brief Deleted copy assignment operator */
    CollisionDetectionModule& operator=(const CollisionDetectionModule&) = delete;

    /** @brief Defaulted move constructor */
    CollisionDetectionModule(CollisionDetectionModule&&) = default;

    /** @brief Defaulted move assignment operator */
    CollisionDetectionModule& operator=(CollisionDetectionModule&&) = default;
    //@}

    /** @name Get methods */
    //@{
    /** @brief Returns raw pointer to the pair list (size = getPairCount()) */
    const IndexPair* getPairList() const;

    /** @brief Returns the current number of active pairs in the neighbor list */
    uint getPairCount() const;

    /** @brief Returns the allocated size of the pair buffers (may exceed getPairCount()) */
    size_t getPairBufferSize() const;

    /** @brief Returns the neighbor list (read-only) */
    const NeighborList<T, M>* getNeighborList() const;
    //@}

    /** @name Methods */
    //@{
    /** @brief Runs the complete collision detection pipeline.
        Pipeline: 1. sortParticles
                  2. updateNeighborList
                  3. computeRelativeTransformations
                  4. filterPairsBV (DEVICE: BV-ON always; BV-OFF only when composites > 0)
                  5. detectCollisionsComponents
                  6. transformContactInfo
        @param rigidBodies   Raw pointer array of RigidBody*
        @param positions     Position buffer (modified in-place by sorter)
        @param orientations  Quaternion buffer (modified in-place by sorter)
        @param velocities    Kinematics buffer (modified in-place by sorter)
        @param torces        Torce buffer (modified in-place by sorter)
        @param bodyTags      BodyTag buffer (encodes shapeId + composite info; modified in-place
                             by sorter; used for filtering)
        @param localPos      Per-component local position offsets (modified in-place by sorter)
        @param localQuat     Per-component local quaternion offsets (modified in-place by sorter)
        @param masterSlot    Per-composite master slot lookup (size = numComposites; rewritten
                             in-place by sortParticles when a sort runs)
        @param contactInfo   ContactInfo world-frame buffer (owned by ComponentManager)
        @param pairList      Pair list buffer (owned by ComponentManager)
        @param counts        Component counts */
    void run(const RigidBody<T>* const*          rigidBodies,
             GrainsMemBuffer<Vector3<T>, M>&     positions,
             GrainsMemBuffer<Quaternion<T>, M>&  orientations,
             GrainsMemBuffer<Kinematics<T>, M>&  velocities,
             GrainsMemBuffer<Torce<T>, M>&       torces,
             GrainsMemBuffer<uint, M>&           bodyTags,
             GrainsMemBuffer<Vector3<T>, M>&     localPos,
             GrainsMemBuffer<Quaternion<T>, M>&  localQuat,
             GrainsMemBuffer<uint, M>&           masterSlot,
             GrainsMemBuffer<ContactInfo<T>, M>& contactInfo,
             GrainsMemBuffer<IndexPair, M>&      pairList,
             ComponentCounts&                    counts);
    //@}

private:
    /** @name Pipeline steps */
    //@{
    /** @brief Resizes per-pair buffers to at least `size` pairs.
        @param pairList    ComponentManager-owned buffer resized in lock-step
        @param size        New minimum capacity
        @param contactInfo ComponentManager-owned buffer resized in lock-step */
    void resizePairBuffers(GrainsMemBuffer<IndexPair, M>&      pairList,
                           GrainsMemBuffer<ContactInfo<T>, M>& contactInfo,
                           uint                                size);

    /** @brief Calls the neighbor list update; resizes buffers and increments the global update
        counter if the list was rebuilt.
        @param positions   Position buffer required by the NL algorithm
        @param pairList    ComponentManager-owned buffer resized if NL grew
        @param contactInfo ComponentManager-owned buffer resized if NL grew
        @param counts      Component counts */
    void updateNeighborList(GrainsMemBuffer<Vector3<T>, M>&     positions,
                            GrainsMemBuffer<IndexPair, M>&      pairList,
                            GrainsMemBuffer<ContactInfo<T>, M>& contactInfo,
                            ComponentCounts&                    counts);

    /** @brief Computes per-pair relative position/quaternion (B expressed in A-local frame).
        @param positions    World position buffer
        @param orientations World quaternion buffer
        @param pairList     Pair list buffer */
    void computeRelativeTransformations(const GrainsMemBuffer<Vector3<T>, M>&    positions,
                                        const GrainsMemBuffer<Quaternion<T>, M>& orientations,
                                        const GrainsMemBuffer<IndexPair, M>&     pairList);

    /** @brief Runs the BV pre-filter.
        Also used as the composite-only filter when BV is OFF but composites exist.
        Writes per-pair pass/fail flags, writes no-contact sentinels for rejected pairs, then
        returns the selected original pair indices from the CUDA backend. The backend owns
        all compaction buffers and the mapped count.
        @param rigidBodies      Raw pointer array of RigidBody*
        @param pairList         Pair list buffer
        @param contactInfoLocal Per-pair local ContactInfo buffer
        @param bodyTags         Raw device pointer to per-component body tags
        @param numComposites    Number of composites (0 = skip composite check) */
    SelectedPairIndices filterPairsBV(const RigidBody<T>* const*           rigidBodies,
                                      const GrainsMemBuffer<IndexPair, M>& pairList,
                                      GrainsMemBuffer<ContactInfo<T>, M>&  contactInfoLocal,
                                      const uint*                          bodyTags,
                                      uint                                 numComposites);

    /** @brief BV pre-filter using world-frame positions/quaternions directly, bypassing the
        computeRelativeTransformations pre-pass. Used when UseRelativeTransformations=false with
        BV enabled.
        @param rigidBodies      Raw pointer array of RigidBody*
        @param positions        World position buffer
        @param orientations     World quaternion buffer
        @param pairList         Pair list buffer
        @param contactInfoLocal Per-pair local ContactInfo buffer
        @param bodyTags         Raw device pointer to per-component body tags
        @param numComposites    Number of composites (0 = skip composite check) */
    SelectedPairIndices filterPairsBV_global(const RigidBody<T>* const*               rigidBodies,
                                             const GrainsMemBuffer<Vector3<T>, M>&    positions,
                                             const GrainsMemBuffer<Quaternion<T>, M>& orientations,
                                             const GrainsMemBuffer<IndexPair, M>&     pairList,
                                             GrainsMemBuffer<ContactInfo<T>, M>& contactInfoLocal,
                                             const uint*                         bodyTags,
                                             uint                                numComposites);

    /** @brief Runs narrow-phase GJK and writes ContactInfo in A-local frame.
        On DEVICE with BV-ON or composites, uses the backend-owned compacted pair-index list.
        On HOST, skips intra-composite pairs inline.
        @param rigidBodies Raw pointer array of RigidBody*
        @param pairList    Pair list buffer
        @param bodyTags    Per-component body tag buffer
        @param counts      Component counts */
    void detectCollisionsComponents(const RigidBody<T>* const*           rigidBodies,
                                    const GrainsMemBuffer<IndexPair, M>& pairList,
                                    const GrainsMemBuffer<uint, M>&      bodyTags,
                                    const ComponentCounts&               counts);

    /** @brief Runs narrow-phase GJK using world-frame positions/quaternions directly, skipping the
        relative-transformation pre-pass.
        On DEVICE with BV-ON or composites, uses the backend-owned compacted pair-index list.
        On HOST, skips intra-composite pairs inline.
        @param rigidBodies  Raw pointer array of RigidBody*
        @param positions    World position buffer
        @param orientations World quaternion buffer
        @param pairList     Pair list buffer
        @param contactInfo  ComponentManager's world-frame ContactInfo buffer
        @param bodyTags     Per-component body tag buffer
        @param counts       Component counts */
    void detectCollisionsComponentsGlobal(const RigidBody<T>* const*               rigidBodies,
                                          const GrainsMemBuffer<Vector3<T>, M>&    positions,
                                          const GrainsMemBuffer<Quaternion<T>, M>& orientations,
                                          const GrainsMemBuffer<IndexPair, M>&     pairList,
                                          GrainsMemBuffer<ContactInfo<T>, M>&      contactInfo,
                                          const GrainsMemBuffer<uint, M>&          bodyTags,
                                          const ComponentCounts&                   counts);

    /** @brief Transforms per-pair ContactInfo from A-local frame to world frame.
        @param positions    World position buffer
        @param orientations World quaternion buffer
        @param pairList     Pair list buffer
        @param contactInfo  ComponentManager's world-frame ContactInfo buffer */
    void transformContactInfo(const GrainsMemBuffer<Vector3<T>, M>&    positions,
                              const GrainsMemBuffer<Quaternion<T>, M>& orientations,
                              const GrainsMemBuffer<IndexPair, M>&     pairList,
                              GrainsMemBuffer<ContactInfo<T>, M>&      contactInfo);

    /** @brief Sorts particles by Morton codes for improved cache efficiency.
        @param positions    Position buffer (reordered in-place)
        @param orientations Quaternion buffer (reordered in-place)
        @param velocities   Kinematics buffer (reordered in-place)
        @param torces       Torce buffer (reordered in-place)
        @param bodyTags     BodyTag buffer (encodes shapeId + composite info; reordered in-place)
        @param localPos     LocalPos buffer (reordered in-place)
        @param localQuat    LocalQuat buffer (reordered in-place)
        @param masterSlot   Master slot buffer (reordered in-place)
        @param counts       Component counts */
    void sortParticles(GrainsMemBuffer<Vector3<T>, M>&    positions,
                       GrainsMemBuffer<Quaternion<T>, M>& orientations,
                       GrainsMemBuffer<Kinematics<T>, M>& velocities,
                       GrainsMemBuffer<Torce<T>, M>&      torces,
                       GrainsMemBuffer<uint, M>&          bodyTags,
                       GrainsMemBuffer<Vector3<T>, M>&    localPos,
                       GrainsMemBuffer<Quaternion<T>, M>& localQuat,
                       GrainsMemBuffer<uint, M>&          masterSlot,
                       const ComponentCounts&             counts);

    /** @brief Builds (or rebuilds) the compact m_shapeData and m_bvData arrays.
        Scans bodyTags to identify unique shapes (by shapeId), builds a representative-slot
        table, then fills one ShapeData and one BVData entry per unique shape.
        Representative slots are built from hostBodyTags on CPU for both backends.
        Backends fill the tables; CUDA privately uploads representative slots.
        @param rigidBodies   Raw pointer array of RigidBody* (slot-indexed)
        @param bodyTagsData  Raw pointer to per-component body tags (host-accessible)
        @param nComponents   Total number of component slots (obstacles + particles) */
    void buildShapeAndBVData(const RigidBody<T>* const* rigidBodies,
                             const uint*                bodyTagsData,
                             uint                       nComponents);
    //@}
};

#endif
