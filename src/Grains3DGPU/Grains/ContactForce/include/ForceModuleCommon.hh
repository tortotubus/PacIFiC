#ifndef _FORCEMODULECOMMON_HH_
#define _FORCEMODULECOMMON_HH_

// Uncomment to enable contact debug output (CPU path only)
// #define GRAINS_CONTACT_DEBUG

#include "Basic.hh"
#include "ContactForceModel.hh"
#include "ContactInfo.hh"
#include "ContactTable.hh"
#include "GrainsParameters.hh"
#include "HookeContactForceModel.hh"
#include "HookeMemoryContactForceModel.hh"
#include "Kinematics.hh"
#include "RigidBody.hh"
#include "Torce.hh"
#include "Vector3.hh"

// =================================================================================================
/** @brief ForceModule common host/device helper functions.

    Header-only inline functions shared by the HOST (CPU loop) and DEVICE (GPU kernel) paths of
    ForceModule.  Functions are templated on the scalar type T and annotated with __HOSTDEVICE__
    or __device__ as appropriate so that a single definition compiles for both paths.

    Functions moved here from ComponentManagerCommon.hh:
      - computeContactForces_common  (single-torce variant, CPU path)
      - computeContactForces_common  (dual-torce variant, GPU intermediate storage path)
      - reduceTorces_common           (__device__ only)
      - addExternalForces_common

    @author A.Yazdani - 2024 - Construction */
// =================================================================================================
/** @name ForceModule common functions between host and device */
//@{
/** @brief Dispatches computeForces with devirtualization when CFType is known at compile time.
    When CFType is HOOKE or HOOKEMEMORY, uses static_cast to call the concrete final method
    directly, eliminating vtable lookup and enabling inlining.  Falls back to virtual dispatch
    when CFType is CF_VIRTUAL (default). */
template <typename T, ContactForceModelType CFType>
__HOSTDEVICE__ static INLINE void dispatchComputeForces(const ContactForceModel<T>* const* CF,
                                                        uint                               cfID,
                                                        const ContactInfo<T>&              ci,
                                                        const Vector3<T>&                  relVel,
                                                        const Vector3<T>&  relAngVel,
                                                        const Vector3<T>&  posA,
                                                        const Vector3<T>&  posB,
                                                        ContactHistory<T>* hist,
                                                        Torce<T>&          tA,
                                                        Torce<T>&          tB)
{
    if constexpr(CFType == HOOKE)
        static_cast<const HookeContactForceModel<T>*>(CF[cfID])
            ->computeForces(ci, relVel, relAngVel, posA, posB, hist, tA, tB);
    else if constexpr(CFType == HOOKEMEMORY)
        static_cast<const HookeMemoryContactForceModel<T>*>(CF[cfID])
            ->computeForces(ci, relVel, relAngVel, posA, posB, hist, tA, tB);
    else
        CF[cfID]->computeForces(ci, relVel, relAngVel, posA, posB, hist, tA, tB);
}

// -------------------------------------------------------------------------------------------------
/** @brief Computes the contact forces and writes directly to per-particle torce array (CPU path).
    @param CF contact force models
    @param pairList list of pairs
    @param contactInfo contact information in the world frame
    @param position position of the components
    @param velocity kinematics of the components
    @param torce per-particle torce array (modified in-place, indexed by component ID)
    @param contactMemory view of contact memory (hash table + history data)
    @param pairID ID of the pair */
template <typename T, ContactForceModelType CFType = CF_VIRTUAL>
__HOSTDEVICE__ static INLINE void
    computeContactForces_common(const ContactForceModel<T>* const* __RESTRICT__ CF,
                                const uint2* __RESTRICT__                       pairList,
                                const ContactInfo<T>* __RESTRICT__              contactInfo,
                                const Vector3<T>* __RESTRICT__                  position,
                                const Kinematics<T>* __RESTRICT__               velocity,
                                Torce<T>* __RESTRICT__                          torce,
                                ContactMemoryView<T>                            contactMemory,
                                const uint                                      pairID)
{
    // Contact Details
    ContactInfo<T> ci = contactInfo[pairID];
    // one load of metadata
    typename ContactInfo<T>::Snapshot snapshot = ci.getSnapshot();
    bool isContact = snapshot.overlapDistance < T(0);  // is in contact / negative distance

    // Compute the forces
    if(isContact)
    {
        const uint2 pair = pairList[pairID];
        const uint  idA  = pair.x;
        const uint  idB  = pair.y;

        // velocities of the components
        const Kinematics<T>& vA(velocity[idA]);
        const Kinematics<T>& vB(velocity[idB]);
        // geometric point of contact
        const Vector3<T>& contactPt(ci.getContactPoint());
        // relative velocity at contact point
        const Vector3<T>& relVel(vA.kinematicsAtPoint(contactPt - position[idA])
                                 - vB.kinematicsAtPoint(contactPt - position[idB]));
        // relative angular velocity
        const Vector3<T>& relAngVel(vA.getAngularComponent() - vB.getAngularComponent());

        // Look up or create contact history entry
        ContactHistory<T>* historyPtr = nullptr;
        if(contactMemory.m_historyData != nullptr)
        {
            uint historyIndex;
            if(contactMemory.findOrInsert(pair, historyIndex))
                historyPtr = &(contactMemory.m_historyData[historyIndex]);
        }

        // note that we will add torce to obstacles as well.
        uint contactForceID = snapshot.contactHash;
        dispatchComputeForces<T, CFType>(CF,
                                         contactForceID,
                                         ci,
                                         relVel,
                                         relAngVel,
                                         position[idA],
                                         position[idB],
                                         historyPtr,
                                         torce[idA],
                                         torce[idB]);
    }
    // reset the distance so we don't compute the torce twice
    ci.setOverlapDistance(T(0));
}

// -------------------------------------------------------------------------------------------------
/** @brief Computes the contact forces (writes to intermediate per-pair storage; GPU path).
    @param CF contact force models
    @param pairList list of pairs
    @param contactInfo contact information in the world frame
    @param position position of the components
    @param velocity kinematics of the components
    @param torceA intermediate torce storage for particle A (indexed by localID)
    @param torceB intermediate torce storage for particle B (indexed by localID)
    @param contactMemory view of contact memory (hash table + history data)
    @param pairID ID of the pair in the full neighbor list (for reading pairList/contactInfo)
    @param localID compact slot index for writing to intermediateTorce buffers
           (equals pairID on the no-compaction path, equals thread index on the compaction path) */
template <typename T, ContactForceModelType CFType = CF_VIRTUAL>
__HOSTDEVICE__ static INLINE void
    computeContactForces_common(const ContactForceModel<T>* const* __RESTRICT__ CF,
                                const uint2* __RESTRICT__                       pairList,
                                const ContactInfo<T>* __RESTRICT__              contactInfo,
                                const Vector3<T>* __RESTRICT__                  position,
                                const Kinematics<T>* __RESTRICT__               velocity,
                                Torce<T>* __RESTRICT__                          torceA,
                                Torce<T>* __RESTRICT__                          torceB,
                                ContactMemoryView<T>                            contactMemory,
                                const uint                                      pairID,
                                const uint                                      localID)
{
    // Contact Details
    ContactInfo<T> ci = contactInfo[pairID];
    // one load of metadata
    typename ContactInfo<T>::Snapshot snapshot = ci.getSnapshot();
    bool isContact = snapshot.overlapDistance < T(0);  // is in contact / negative distance

    // Compute the forces
    if(isContact)  // On device path, this check is redundant.
    {
        const uint2 pair = pairList[pairID];
        const uint  idA  = pair.x;
        const uint  idB  = pair.y;

        // velocities of the components
        const Kinematics<T>& vA(velocity[idA]);
        const Kinematics<T>& vB(velocity[idB]);
        // geometric point of contact
        const Vector3<T>& contactPt(ci.getContactPoint());
        // relative velocity at contact point
        const Vector3<T>& relVel(vA.kinematicsAtPoint(contactPt - position[idA])
                                 - vB.kinematicsAtPoint(contactPt - position[idB]));
        // relative angular velocity
        const Vector3<T>& relAngVel(vA.getAngularComponent() - vB.getAngularComponent());

        // Look up or create contact history entry
        ContactHistory<T>* historyPtr = nullptr;
        if(contactMemory.m_historyData != nullptr)
        {
            uint historyIndex;
            if(contactMemory.findOrInsert(pair, historyIndex))
                historyPtr = &(contactMemory.m_historyData[historyIndex]);
        }

        // Reset intermediate torces before accumulating -- indexed by localID (compact slot).
        torceA[localID].reset();
        torceB[localID].reset();

        // note that we will add torce to obstacles as well.
        uint contactForceID = snapshot.contactHash;
        dispatchComputeForces<T, CFType>(CF,
                                         contactForceID,
                                         ci,
                                         relVel,
                                         relAngVel,
                                         position[idA],
                                         position[idB],
                                         historyPtr,
                                         torceA[localID],
                                         torceB[localID]);
    }
    // reset the distance so we don't compute the torce twice
    ci.setOverlapDistance(T(0));
}

// -------------------------------------------------------------------------------------------------
/** @brief Reduces per-pair intermediate torces to per-particle torces (DEVICE kernel helper).
    @param pairList list of pairs
    @param intermediateTorceA intermediate torce storage for particle A
    @param intermediateTorceB intermediate torce storage for particle B
    @param torce final per-particle torce array (accumulated atomically)
    @param pairID ID of the pair in the full neighbor list (for reading pairList)
    @param localID compact slot index for reading from intermediateTorce buffers */
template <typename T>
__device__ static INLINE void reduceTorces_common(const uint2* __RESTRICT__    pairList,
                                                  const Torce<T>* __RESTRICT__ intermediateTorceA,
                                                  const Torce<T>* __RESTRICT__ intermediateTorceB,
                                                  Torce<T>* __RESTRICT__       torce,
                                                  const uint                   pairID,
                                                  const uint                   localID)
{
    const uint2     pair = pairList[pairID];
    const uint      idA  = pair.x;
    const uint      idB  = pair.y;
    const Torce<T>& tA   = intermediateTorceA[localID];
    const Torce<T>& tB   = intermediateTorceB[localID];

    // Atomically accumulate all components (6 atomics per particle)
    torce[idA].addTorceAtomic(tA);
    torce[idB].addTorceAtomic(tB);
}

// -------------------------------------------------------------------------------------------------
/** @brief Adds gravity to the component.
    @param g the gravitational acceleration vector
    @param rigidBody the rigid body of the component
    @param torce the torce acting on the component
    @param cID the ID of the component */
template <typename T>
__HOSTDEVICE__ static INLINE void addExternalForces_common(const Vector3<T>&          g,
                                                           const RigidBody<T>* const* rigidBody,
                                                           Torce<T>*                  torce,
                                                           const uint                 cID)
{
    const RigidBody<T>* rb   = rigidBody[cID];
    const T             mass = rb->getMass();
    // Adding the gravitational force to the torce
    torce[cID].addForce(mass * g);
}
//@}

#endif
