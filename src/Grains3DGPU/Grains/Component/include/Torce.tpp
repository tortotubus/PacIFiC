#ifndef TORCE_TPP
#define TORCE_TPP

template <typename T>
__HOSTDEVICE__ Torce<T>::Torce()
    : m_torque(0, 0, 0)
    , m_force(0, 0, 0)
{
}

template <typename T>
__HOSTDEVICE__ Torce<T>::Torce(const Vector3<T>& torque, const Vector3<T>& force)
    : m_torque(torque)
    , m_force(force)
{
}

template <typename T>
__HOSTDEVICE__ Torce<T>::~Torce()
{
}

template <typename T>
__HOSTDEVICE__ Vector3<T> Torce<T>::getTorque() const
{
    return m_torque;
}

template <typename T>
__HOSTDEVICE__ Vector3<T> Torce<T>::getForce() const
{
    return m_force;
}

template <typename T>
__HOSTDEVICE__ void Torce<T>::setTorque(const Vector3<T>& torque)
{
    m_torque = torque;
}

template <typename T>
__HOSTDEVICE__ void Torce<T>::setForce(const Vector3<T>& force)
{
    m_force = force;
}

template <typename T>
__HOSTDEVICE__ void Torce<T>::reset()
{
    m_torque.reset();
    m_force.reset();
}

template <typename T>
__HOSTDEVICE__ void Torce<T>::addTorque(const Vector3<T>& torque)
{
    for(uint i = 0; i < 3; ++i)
        m_torque[i] += torque[i];
}

template <typename T>
__HOSTDEVICE__ void Torce<T>::addForce(const Vector3<T>& force)
{
    for(uint i = 0; i < 3; ++i)
        m_force[i] += force[i];
}

template <typename T>
__HOSTDEVICE__ void Torce<T>::addForce(const Vector3<T>& force, const Vector3<T>& point)
{
    addForce(force);
    m_torque[X] += point[Y] * force[Z] - point[Z] * force[Y];
    m_torque[Y] += point[Z] * force[X] - point[X] * force[Z];
    m_torque[Z] += point[X] * force[Y] - point[Y] * force[X];
}

#endif
