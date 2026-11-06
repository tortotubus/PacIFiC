#ifndef VECTOR3_TPP
#define VECTOR3_TPP

template <typename T>
__HOSTDEVICE__ Vector3<T>::Vector3(T def) noexcept
{
    m_comp[X] = m_comp[Y] = m_comp[Z] = def;
}

template <typename T>
__HOSTDEVICE__ Vector3<T>::Vector3(const T* buffer) noexcept
{
    setValue(buffer);
}

template <typename T>
__HOSTDEVICE__ Vector3<T>::Vector3(T x, T y, T z) noexcept
{
    setValue(x, y, z);
}

template <typename T>
__HOSTDEVICE__ Vector3<T>::Vector3(const Vector3<T>& vector) noexcept
{
    setValue(vector.m_comp);
}

template <typename T>
__HOSTDEVICE__ Vector3<T>& Vector3<T>::operator=(const Vector3<T>& vector) noexcept
{
    if(&vector != this)
        setValue(vector.m_comp);
    return *this;
}

template <typename T>
__HOSTDEVICE__ Vector3<T>::Vector3(Vector3<T>&& vector) noexcept
{
    setValue(vector.m_comp);
    vector.reset();
}

template <typename T>
__HOSTDEVICE__ Vector3<T>& Vector3<T>::operator=(Vector3<T>&& vector) noexcept
{
    if(&vector != this)
    {
        setValue(vector.m_comp);
        vector.reset();
    }
    return *this;
}

template <typename T>
__HOSTDEVICE__ Vector3<T>::~Vector3() noexcept
{
}

template <typename T>
__HOSTDEVICE__ const T* Vector3<T>::getBuffer() const noexcept
{
    return m_comp;
}

template <typename T>
__HOSTDEVICE__ void Vector3<T>::setValue(const T* buffer) noexcept
{
    m_comp[X] = buffer[X];
    m_comp[Y] = buffer[Y];
    m_comp[Z] = buffer[Z];
}

template <typename T>
__HOSTDEVICE__ void Vector3<T>::setValue(const T x, const T y, const T z) noexcept
{
    m_comp[X] = x;
    m_comp[Y] = y;
    m_comp[Z] = z;
}

template <typename T>
__HOSTDEVICE__ void Vector3<T>::normalize() noexcept
{
    const T magnitude
        = sqrt(m_comp[X] * m_comp[X] + m_comp[Y] * m_comp[Y] + m_comp[Z] * m_comp[Z]);
    m_comp[X] /= magnitude;
    m_comp[Y] /= magnitude;
    m_comp[Z] /= magnitude;
}

template <typename T>
__HOSTDEVICE__ Vector3<T> Vector3<T>::normalized() const noexcept
{
    const T magnitude
        = sqrt(m_comp[X] * m_comp[X] + m_comp[Y] * m_comp[Y] + m_comp[Z] * m_comp[Z]);
    if(magnitude > EPS<T>)
        return Vector3<T>(m_comp[X] / magnitude, m_comp[Y] / magnitude, m_comp[Z] / magnitude);
    return Vector3<T>();
}

template <typename T>
__HOSTDEVICE__ void Vector3<T>::reset() noexcept
{
    m_comp[X] = m_comp[Y] = m_comp[Z] = T(0);
}

template <typename T>
__HOSTDEVICE__ const T& Vector3<T>::operator[](size_t i) const noexcept
{
    return m_comp[i];
}

template <typename T>
__HOSTDEVICE__ T& Vector3<T>::operator[](size_t i) noexcept
{
    return m_comp[i];
}

#endif
