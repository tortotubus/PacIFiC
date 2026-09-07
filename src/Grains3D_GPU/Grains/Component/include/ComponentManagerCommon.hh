#ifndef _COMPONENTMANAGERCOMMON_HH_
#define _COMPONENTMANAGERCOMMON_HH_

#include "Basic.hh"
#include "GrainsParameters.hh"
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
/** @brief ComponentManager common functions between host and device.

    This is a header-only file that contains common functions between CPU/GPU
    for the ComponentManager class. The functions are templated to allow
    for flexibility in usage. The functions are marked as inline to allow for
    better optimization by the compiler.

    Force computation helpers (computeContactForces_common, reduceTorces_common,
    addExternalForces_common) have been moved to ForceModuleCommon.hh.

    @author A.Yazdani - 2025 - Construction */
// =================================================================================================
/** @name ComponentManager common functions between host and device */
//@{
/** @brief Applies one prescribed obstacle-motion event.
    Returns true only when the event is active at currentTime and produces motion in this step.
    The translational update uses Vector3::atomicAdd, which falls back to regular addition on
    host and uses CUDA atomics on device.
    @param position      array of obstacle positions
    @param quaternion    array of obstacle quaternions
    @param ev            obstacle motion event
    @param currentTime   current simulation time
    @param dt            time step */
template <typename T>
__HOSTDEVICE__ static INLINE bool moveObstacle_common(Vector3<T>*                   position,
                                                      Quaternion<T>*                quaternion,
                                                      Kinematics<T>*                velocity,
                                                      const ObstacleMotionEvent<T>& ev,
                                                      const T                       currentTime,
                                                      const T                       dt)
{
    if(currentTime < ev.tStart || currentTime >= ev.tEnd)
        return false;

    const T linSq = norm2(ev.linearVelocity);
    const T angSq = norm2(ev.angularVelocity);
    if(linSq < EPS<T> * EPS<T> && angSq < EPS<T> * EPS<T>)
        return false;

    const uint id = ev.obstacleId;
    velocity[id].setTranslationalComponent(ev.linearVelocity);
    velocity[id].setAngularComponent(ev.angularVelocity);
    // on host, it falls back to regular addition.
    position[id].atomicAdd(dt * ev.linearVelocity);

    const T nOmega = sqrt(angSq);
    if(nOmega > EPS<T>)
    {
        const T             half_angle = nOmega * dt / T(2);
        const T             s          = sin(half_angle) / nOmega;
        const Quaternion<T> dq(s * ev.angularVelocity, cos(half_angle));
        quaternion[id] = dq * quaternion[id];
        const T qn     = norm(quaternion[id]);
        if(qn > EPS<T>)
            quaternion[id] *= (T(1) / qn);
    }

    return true;
}

// -------------------------------------------------------------------------------------------------
/** @brief Moves a component using the given time integration method.
    @param TI                  the time integrator
    @param rigidBody           the rigid body of the components
    @param position            world-frame positions (standalone updated in place)
    @param compositePosition   composite CM world position (one per composite; updated for masters)
    @param quaternion          world-frame quaternions (standalone updated in place)
    @param compositeQuaternion composite CM world quaternion (one per composite; updated for
   masters)
    @param velocity            per-component kinematics (standalone updated; composites via
   updateSubBody)
    @param compositeVelocity composite CM velocity (one per composite; updated for masters via Move)
    @param torce               per-component torce (standalone reset; composite entries unused)
    @param compositeTorce      composite net force/torque (one per composite; reset for masters)
    @param bodyTag             per-component body tag
    @param cID                 the ID of the component */
template <typename T>
__HOSTDEVICE__ static INLINE void moveParticles_common(const TimeIntegrator<T>* const* TI,
                                                       const RigidBody<T>* const*      rigidBody,
                                                       Vector3<T>*                     position,
                                                       Vector3<T>*    compositePosition,
                                                       Quaternion<T>* quaternion,
                                                       Quaternion<T>* compositeQuaternion,
                                                       Kinematics<T>* velocity,
                                                       Kinematics<T>* compositeVelocity,
                                                       Torce<T>*      torce,
                                                       Torce<T>*      compositeTorce,
                                                       const uint*    bodyTag,
                                                       const uint     cID,
                                                       const T        maxDisplacementSq)
{
    // Rigid body
    const RigidBody<T>* rb = rigidBody[cID];

    const uint tag          = bodyTag[cID];
    const bool isCompMaster = isSubBody(tag) && (getSubBodyLocalIdx(tag) == 0u);
    const uint cIdx         = isCompMaster ? getCompositeIdx(tag) : 0u;

    // Select the orientation used for the Euler equations in the principal-axis frame
    const Quaternion<T>& Q = isCompMaster ? compositeQuaternion[cIdx] : quaternion[cID];

    // Select current kinematics and accumulated torce (composite or per-body)
    Kinematics<T>& Kin   = isCompMaster ? compositeVelocity[cIdx] : velocity[cID];
    Torce<T>&      Torce = isCompMaster ? compositeTorce[cIdx] : torce[cID];

    // Computing momentums in the space-fixed coordinate
    const Kinematics<T>& momentum = rb->computeMomentum(Kin.getAngularComponent(), Torce, Q);
    // Reset the torce that was just consumed
    Torce.reset();
    // Move: updates Kin in-place and returns position/rotation increments
    Vector3<T>    transMotion;
    Quaternion<T> rotMotion;
    TI[0]->Move(momentum, Kin, transMotion, rotMotion);

    // Stability check: abort if any particle displaces more than one cell width per step.
    // A particle jumping over an entire cell makes the neighbor list blind to new contacts.
    GAssert(norm2(transMotion) <= maxDisplacementSq,
            "STABILITY FAILED: particle",
            cID,
            "displacement",
            norm(transMotion),
            "exceeds cell width. Simulation is numerically unstable.");

    if(isCompMaster)
    {
        compositePosition[cIdx] += transMotion;
        compositeQuaternion[cIdx] = rotMotion * compositeQuaternion[cIdx];
        T qn                      = norm(compositeQuaternion[cIdx]);
        if(qn > EPS<T>)
            compositeQuaternion[cIdx] *= (T(1) / qn);
    }
    else
    {
        position[cID] += transMotion;
        quaternion[cID] = rotMotion * quaternion[cID];
        T qn            = norm(quaternion[cID]);
        if(qn > EPS<T>)
            quaternion[cID] *= (T(1) / qn);
    }
}

// -------------------------------------------------------------------------------------------------
/** @brief Performs the second velocity half-kick for split-step schemes (KDK leapfrog Step 3).
    @param TI                  the time integrator
    @param rigidBody           the rigid body of the components
    @param quaternion          world-frame quaternions (standalone; also used if not a master)
    @param compositeQuaternion composite CM quaternion (used for Euler equations of masters)
    @param velocity            per-component kinematics (standalone updated; composites via
   updateSubBody)
    @param compositeVelocity composite CM velocity (updated for masters)
    @param torce               per-component torce (standalone read, NOT reset)
    @param compositeTorce      composite net force/torque (read but NOT reset)
    @param bodyTag             per-component body tag
    @param cID                 the ID of the component */
template <typename T>
__HOSTDEVICE__ static INLINE void advanceVelocity_common(const TimeIntegrator<T>* const* TI,
                                                         const RigidBody<T>* const*      rigidBody,
                                                         const Quaternion<T>*            quaternion,
                                                         const Quaternion<T>* compositeQuaternion,
                                                         Kinematics<T>*       velocity,
                                                         Kinematics<T>*       compositeVelocity,
                                                         const Torce<T>*      torce,
                                                         const Torce<T>*      compositeTorce,
                                                         const uint*          bodyTag,
                                                         const uint           cID)
{
    const RigidBody<T>* rb = rigidBody[cID];

    const uint tag          = bodyTag[cID];
    const bool isCompMaster = isSubBody(tag) && (getSubBodyLocalIdx(tag) == 0u);
    const uint cIdx         = isCompMaster ? getCompositeIdx(tag) : 0u;

    const Quaternion<T>& Q            = isCompMaster ? compositeQuaternion[cIdx] : quaternion[cID];
    Kinematics<T>&       Kin          = isCompMaster ? compositeVelocity[cIdx] : velocity[cID];
    const Torce<T>&      Torce        = isCompMaster ? compositeTorce[cIdx] : torce[cID];
    const Kinematics<T>  acceleration = rb->computeMomentum(Kin.getAngularComponent(), Torce, Q);
    TI[0]->AdvanceVelocity(acceleration, Kin);
}
//@}

#endif
