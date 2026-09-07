// TODO: CHANGE THE FORMAT FROM HH TO CUH LATER.
#ifndef _COMPONENTMANAGERGPU_KERNLES_CUH_
#define _COMPONENTMANAGERGPU_KERNLES_CUH_

#include <cuda_runtime.h>

#include "ComponentManager.hh"
#include "ComponentManagerCommon.hh"
#include "Kinematics.hh"
#include "ObstacleLoading.hh"
#include "Quaternion.hh"
#include "QuaternionMath.hh"
#include "RigidBody.hh"
#include "TimeIntegrator.hh"
#include "Torce.hh"
#include "Vector3.hh"
#include "VectorMath.hh"

// =================================================================================================
/** @brief GPU kernels for the ComponentManagerGPU class (particle motion).

    Contact-force kernels have been moved to ForceModule_Kernels.hh.

    @author A.Yazdani - 2024 - Construction */
// =================================================================================================
/** @name ComponentManagerGPU_Kernels : External methods */
//@{
/** @brief Integrates obstacle positions and orientations by one time step.
    One thread per ObstacleMotionEvent entry.  The thread checks whether @p currentTime
    falls within its event's time window and, if so, integrates the indexed obstacle via
    moveObstacle_common(). Obstacles are updated through atomicAdd so the kernel is
    formally data-race-free even if the non-overlap guarantee were ever violated.
    @param events       per-event records (one per obstacle x interval pair)
    @param numEvents    total number of events
    @param currentTime  simulation time at the start of this step
    @param dt           time step length
    @param positions    world-frame positions  (in/out)
    @param quaternions  world-frame quaternions (in/out)
    @param movedFlag    pinned host-mapped integer to indicate if any obstacle moved */
template <typename T>
__GLOBAL__ void moveObstacles_Kernel(const ObstacleMotionEvent<T>* events,
                                     uint                          numEvents,
                                     T                             currentTime,
                                     T                             dt,
                                     Vector3<T>*                   positions,
                                     Quaternion<T>*                quaternions,
                                     Kinematics<T>*                velocities,
                                     int*                          movedFlag)
{
    const uint i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i >= numEvents)
        return;

    if(moveObstacle_common(positions, quaternions, velocities, events[i], currentTime, dt))
        atomicOr(movedFlag, 1);
}

// -------------------------------------------------------------------------------------------------
/** @brief Updates the position and velocities of particles (KDK Leapfrog Step 1).
    For composite masters: advances compositePosition/compositeQuaternion/compositeVelocity
    using compositeTorce; positions/velocities of all sub-bodies are then re-derived by
    updateSubBodies_Kernel.  For standalone particles: advances position/quaternion/velocity.
    Non-master sub-bodies are skipped entirely.
    @param TI                  time integrator scheme
    @param rigidBody           array of rigid bodies for components
    @param position            world-frame positions
    @param compositePosition   composite CM world position (one per composite)
    @param quaternion          array of components quaternions
    @param compositeQuaternion composite CM world quaternion (one per composite)
    @param velocity            array of components velocities
    @param compositeVelocity composite CM velocity (one per composite; updated in place)
    @param torce               array of components torces (standalone; composite entries unused)
    @param compositeTorce      composite net force/torque (one per composite; reset after use)
    @param bodyTag             per-component body tag
    @param nObstacles          number of obstacles
    @param nParticles          number of particles */
template <typename T>
__GLOBAL__ void moveParticles_Kernel(const TimeIntegrator<T>* const* TI,
                                     const RigidBody<T>* const*      rigidBody,
                                     Vector3<T>*                     position,
                                     Vector3<T>*                     compositePosition,
                                     Quaternion<T>*                  quaternion,
                                     Quaternion<T>*                  compositeQuaternion,
                                     Kinematics<T>*                  velocity,
                                     Kinematics<T>*                  compositeVelocity,
                                     Torce<T>*                       torce,
                                     Torce<T>*                       compositeTorce,
                                     const uint*                     bodyTag,
                                     const uint                      nObstacles,
                                     const uint                      nParticles,
                                     const T                         maxDisplacementSq)
{
    uint pID = blockIdx.x * blockDim.x + threadIdx.x;

    if(pID >= nParticles)
        return;

    const uint cID = nObstacles + pID;
    // Skip non-master sub-bodies -- their transforms and velocities are set by
    // updateSubBodies_Kernel
    if(isSubBody(bodyTag[cID]) && getSubBodyLocalIdx(bodyTag[cID]) != 0u)
        return;

    moveParticles_common(TI,
                         rigidBody,
                         position,
                         compositePosition,
                         quaternion,
                         compositeQuaternion,
                         velocity,
                         compositeVelocity,
                         torce,
                         compositeTorce,
                         bodyTag,
                         cID,
                         maxDisplacementSq);
}

// -------------------------------------------------------------------------------------------------
/** @brief Performs the second velocity half-kick for split-step schemes (KDK Leapfrog Step 3).
    For composite masters: updates compositeVelocity using compositeTorce (NOT reset).
    For standalone particles: updates velocity using torce (NOT reset).
    Sub-body velocities are re-derived by updateSubBodyPositions_Kernel after this call.
    @param TI                  time integrator scheme
    @param rigidBody           array of rigid bodies for components
    @param quaternion          world-frame quaternions (standalone)
    @param compositeQuaternion composite CM quaternion (one per composite)
    @param velocity            per-component kinematics (standalone updated)
    @param compositeVelocity   composite CM velocity (one per composite; updated in place)
    @param torce               per-component torce (standalone read, NOT reset)
    @param compositeTorce      composite net force/torque (one per composite; read, NOT reset)
    @param bodyTag             per-component body tag
    @param nObstacles          number of obstacles
    @param nParticles          number of particles */
template <typename T>
__GLOBAL__ void advanceVelocity_Kernel(const TimeIntegrator<T>* const* TI,
                                       const RigidBody<T>* const*      rigidBody,
                                       const Quaternion<T>*            quaternion,
                                       const Quaternion<T>*            compositeQuaternion,
                                       Kinematics<T>*                  velocity,
                                       Kinematics<T>*                  compositeVelocity,
                                       const Torce<T>*                 torce,
                                       const Torce<T>*                 compositeTorce,
                                       const uint*                     bodyTag,
                                       const uint                      nObstacles,
                                       const uint                      nParticles)
{
    uint pID = blockIdx.x * blockDim.x + threadIdx.x;

    if(pID >= nParticles)
        return;

    const uint cID = nObstacles + pID;
    // Skip non-master sub-bodies -- their velocities are re-derived by updateSubBodies_Kernel
    if(isSubBody(bodyTag[cID]) && getSubBodyLocalIdx(bodyTag[cID]) != 0u)
        return;

    advanceVelocity_common(TI,
                           rigidBody,
                           quaternion,
                           compositeQuaternion,
                           velocity,
                           compositeVelocity,
                           torce,
                           compositeTorce,
                           bodyTag,
                           cID);
}
// -------------------------------------------------------------------------------------------------
/** @brief Slaves all sub-body world transforms (position/quaternion) and velocities to their
    composite CM state (compositePosition / compositeQuaternion / compositeVelocity).
    @param localPos            per-component local position offset in composite body frame
    @param localQuat           per-component local quaternion offset in composite body frame
    @param position            world-frame positions (write all sub-bodies)
    @param compositePosition   CM world position per composite (indexed by compositeIdx)
    @param quaternion          world-frame quaternions (write all sub-bodies)
    @param compositeQuaternion CM world quaternion per composite (indexed by compositeIdx)
    @param velocity            per-component kinematics (write v_k = v_CM + omega x r_k)
    @param compositeVelocity   CM velocity per composite (indexed by compositeIdx)
    @param bodyTag             per-component body tag
    @param nComponents         total number of components (obstacles + particles) */
template <typename T>
__GLOBAL__ void updateSubBodies_Kernel(const Vector3<T>*    localPos,
                                       const Quaternion<T>* localQuat,
                                       Vector3<T>*          position,
                                       const Vector3<T>*    compositePosition,
                                       Quaternion<T>*       quaternion,
                                       const Quaternion<T>* compositeQuaternion,
                                       Kinematics<T>*       velocity,
                                       const Kinematics<T>* compositeVelocity,
                                       const uint*          bodyTag,
                                       const uint           nComponents)
{
    uint cID = blockIdx.x * blockDim.x + threadIdx.x;
    if(cID >= nComponents)
        return;

    const uint tag = bodyTag[cID];
    if(!isSubBody(tag))
        return;

    const uint cIdx = getCompositeIdx(tag);
    // world_pos  = compositePosition + rotate(compositeQuaternion, localPos)
    // world_quat = compositeQuaternion * localQuat
    position[cID]   = compositePosition[cIdx] + (compositeQuaternion[cIdx] >> localPos[cID]);
    quaternion[cID] = compositeQuaternion[cIdx] * localQuat[cID];
    T qn            = norm(quaternion[cID]);
    if(qn > EPS<T>)
        quaternion[cID] *= (T(1) / qn);
    // velocity of sub-body CM: v_k = v_CM + omega x r_k
    const Vector3<T>     r     = position[cID] - compositePosition[cIdx];
    const Kinematics<T>& compK = compositeVelocity[cIdx];
    velocity[cID]
        = Kinematics<T>(compK.getTranslationalComponent() + (compK.getAngularComponent() ^ r),
                        compK.getAngularComponent());
}
//@}

#endif
