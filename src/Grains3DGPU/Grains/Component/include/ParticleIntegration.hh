#ifndef GRAINS_PARTICLE_INTEGRATION_HH
#define GRAINS_PARTICLE_INTEGRATION_HH

#include "ComponentManagerCommon.hh"

// Non-owning arrays in a single memory space. Obstacles occupy [0, begin); the
// particle range is [begin, end). Each index can be advanced independently.
template <typename T>
struct ParticleVelocityView
{
    const TimeIntegrator<T>* const* integrator;
    const RigidBody<T>* const*      rigidBody;
    const Quaternion<T>*            quaternion;
    Kinematics<T>*                  velocity;
    const Torce<T>*                 torce;
    const uint*                     bodyTag;
    uint                            begin;
    uint                            end;

    __HOSTDEVICE__ bool integrates(uint i) const
    {
        return i >= begin && i < end
               && (!isSubBody(bodyTag[i]) || getSubBodyLocalIdx(bodyTag[i]) == 0u);
    }

    __HOSTDEVICE__ void advanceVelocity(uint i) const
    {
        if(integrates(i))
            advanceVelocity_common(integrator, rigidBody, quaternion, velocity, torce, i);
    }
};

template <typename T>
struct ParticleIntegrationView
{
    const TimeIntegrator<T>* const* integrator;
    const RigidBody<T>* const*      rigidBody;
    Vector3<T>*                     position;
    Quaternion<T>*                  quaternion;
    Kinematics<T>*                  velocity;
    Torce<T>*                       torce;
    const uint*                     bodyTag;
    uint                            begin;
    uint                            end;

    __HOSTDEVICE__ ParticleVelocityView<T> velocityView() const
    {
        return {integrator, rigidBody, quaternion, velocity, torce, bodyTag, begin, end};
    }

    __HOSTDEVICE__ bool integrates(uint i) const
    {
        return velocityView().integrates(i);
    }

    __HOSTDEVICE__ void move(uint i) const
    {
        if(integrates(i))
            moveParticles_common(integrator, rigidBody, position, quaternion, velocity, torce, i);
    }
};

// Non-owning composite transforms. Backends choose the iteration range and
// normalization tolerance to preserve their existing behavior.
template <typename T>
struct SubBodyTransformView
{
    Vector3<T>*          position;
    Quaternion<T>*       quaternion;
    const Vector3<T>*    localPos;
    const Quaternion<T>* localQuat;
    const uint*          masterSlot;
    const uint*          bodyTag;
    uint                 begin;
    uint                 end;

    __HOSTDEVICE__ void update(uint i, T normalizationTolerance) const
    {
        if(i < begin || i >= end)
            return;
        const uint tag = bodyTag[i];
        if(!isSubBody(tag) || getSubBodyLocalIdx(tag) == 0u)
            return;
        const uint master = masterSlot[getCompositeIdx(tag)];
        position[i]       = position[master] + (quaternion[master] >> localPos[i]);
        quaternion[i]     = quaternion[master] * localQuat[i];
        const T qn        = norm(quaternion[i]);
        if(qn > normalizationTolerance)
            quaternion[i] *= T(1) / qn;
    }
};

#endif
