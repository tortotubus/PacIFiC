#ifndef TORCE_DEVICE_OPS_CUH
#define TORCE_DEVICE_OPS_CUH

#include "Torce.hh"

namespace grains::cuda
{
struct Vector3Ops
{
    template <typename T>
    __DEVICE__ static void atomicAdd(Vector3<T>& target, const Vector3<T>& value) noexcept
    {
        ::atomicAdd(target.m_comp + 0, value.m_comp[0]);
        ::atomicAdd(target.m_comp + 1, value.m_comp[1]);
        ::atomicAdd(target.m_comp + 2, value.m_comp[2]);
    }
};

struct TorceOps
{
    template <typename T>
    __DEVICE__ static void addTorqueAtomic(Torce<T>& target, const Vector3<T>& torque)
    {
        Vector3Ops::atomicAdd(target.m_torque, torque);
    }

    template <typename T>
    __DEVICE__ static void addForceAtomic(Torce<T>& target, const Vector3<T>& force)
    {
        Vector3Ops::atomicAdd(target.m_force, force);
    }

    template <typename T>
    __DEVICE__ static void addTorceAtomic(Torce<T>& target, const Torce<T>& value)
    {
        Vector3Ops::atomicAdd(target.m_force, value.m_force);
        Vector3Ops::atomicAdd(target.m_torque, value.m_torque);
    }
};
}  // namespace grains::cuda

#endif
