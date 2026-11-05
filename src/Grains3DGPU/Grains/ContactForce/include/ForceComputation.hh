#ifndef GRAINS_FORCE_COMPUTATION_HH
#define GRAINS_FORCE_COMPUTATION_HH

#include "BodyTag.hh"
#include "ForceModuleCommon.hh"

// Non-owning views of arrays in one memory space. Ownership and scheduling stay
// with the backend; contact laws and composite moment arms are shared.
template <typename T>
struct ContactForceView
{
    const ContactForceModel<T>* const* models;
    const IndexPair*                   pairs;
    const ContactInfo<T>*              contacts;
    const Vector3<T>*                  positions;
    const Kinematics<T>*               velocities;
    Torce<T>*                          torces;
    ContactMemoryView<T>               memory;
    uint                               count;

    __HOSTDEVICE__ void compute(uint i) const
    {
        if(i < count)
            computeContactForces_common(models,
                                        pairs,
                                        contacts,
                                        positions,
                                        velocities,
                                        torces,
                                        memory,
                                        i);
    }
};

template <typename T>
struct PairForceView
{
    ContactForceView<T> contact;
    Torce<T>*           first;
    Torce<T>*           second;

    __HOSTDEVICE__ void compute(uint i) const
    {
        if(i >= contact.count)
            return;
        // Inactive pairs must not carry forces from the previous step into reduction.
        first[i].reset();
        second[i].reset();
        computeContactForces_common(contact.models,
                                    contact.pairs,
                                    contact.contacts,
                                    contact.positions,
                                    contact.velocities,
                                    first,
                                    second,
                                    contact.memory,
                                    i);
    }
};

template <typename T>
struct ExternalForceView
{
    const RigidBody<T>* const* rigidBodies;
    Torce<T>*                  torces;
    T                          gx, gy, gz;
    uint                       begin, end;

    __HOSTDEVICE__ void add(uint i) const
    {
        if(i >= begin && i < end)
            addExternalForces_common(Vector3<T>(gx, gy, gz), rigidBodies, torces, i);
    }
};

template <typename T>
struct CompositeForceView
{
    Torce<T>*         torces;
    const Vector3<T>* positions;
    const uint*       masterSlots;
    const uint*       bodyTags;
    uint              begin, end;

    __HOSTDEVICE__ bool
        contribution(uint i, uint& master, Vector3<T>& force, Vector3<T>& torque) const
    {
        if(i < begin || i >= end)
            return false;
        const uint tag = bodyTags[i];
        if(!isSubBody(tag) || getSubBodyLocalIdx(tag) == 0u)
            return false;
        master               = masterSlots[getCompositeIdx(tag)];
        const Vector3<T> arm = positions[i] - positions[master];
        force                = torces[i].getForce();
        torque               = torces[i].getTorque() + (arm ^ force);
        return true;
    }
};

#endif
