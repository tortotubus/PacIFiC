#ifndef GRAINS_GEOMETRY_BOX_IMPL_HH
#define GRAINS_GEOMETRY_BOX_IMPL_HH

#include "Box.hh"
#include "VectorMath.hh"

// -------------------------------------------------------------------------------------------------
// Constructor with half edge length as input parameters
template <typename T>
__HOSTDEVICE__ Box<T>::Box(T x, T y, T z)
    : m_extent(Vector3<T>(x, y, z))
{
}
// -------------------------------------------------------------------------------------------------
// Constructor with a vector containing the edge half-lengths
template <typename T>
__HOSTDEVICE__ Box<T>::Box(const Vector3<T>& extent_)
    : m_extent(extent_)
{
}
// -------------------------------------------------------------------------------------------------
// Destructor
template <typename T>
__HOSTDEVICE__ Box<T>::~Box()
{
}
// -------------------------------------------------------------------------------------------------
// Returns the convex type
template <typename T>
__HOSTDEVICE__ ConvexType Box<T>::getConvexType() const
{
    return (BOX);
}
// -------------------------------------------------------------------------------------------------
// Gets values of the edge length
template <typename T>
__HOSTDEVICE__ Vector3<T> Box<T>::getExtent() const
{
    return (m_extent);
}
// -------------------------------------------------------------------------------------------------
// Sets values of the edge length
template <typename T>
__HOSTDEVICE__ void Box<T>::setExtent(T x, T y, T z)
{
    m_extent = Vector3<T>(x, y, z);
}
// -------------------------------------------------------------------------------------------------
// Returns a clone of the box
template <typename T>
__HOSTDEVICE__ Convex<T>* Box<T>::clone() const
{
    return (new Box<T>(m_extent));
}
// -------------------------------------------------------------------------------------------------
// Returns the volume of the box
template <typename T>
__HOSTDEVICE__ T Box<T>::computeVolume() const
{
    return (T(8) * m_extent[X] * m_extent[Y] * m_extent[Z]);
}
// -------------------------------------------------------------------------------------------------
// Computes the diagonal inertia tensor
template <typename T>
__HOSTDEVICE__ void Box<T>::computeInertia(T (&inertia)[3]) const
{
    // Diagonal components: Ixx, Iyy, Izz
    inertia[0] = T(8) * m_extent[X] * m_extent[Y] * m_extent[Z]
                 * (m_extent[Y] * m_extent[Y] + m_extent[Z] * m_extent[Z]) / T(3);
    inertia[1] = T(8) * m_extent[X] * m_extent[Y] * m_extent[Z]
                 * (m_extent[X] * m_extent[X] + m_extent[Z] * m_extent[Z]) / T(3);
    inertia[2] = T(8) * m_extent[X] * m_extent[Y] * m_extent[Z]
                 * (m_extent[Y] * m_extent[Y] + m_extent[X] * m_extent[X]) / T(3);
}
// -------------------------------------------------------------------------------------------------
// Returns the circumscribed radius of the box
template <typename T>
__HOSTDEVICE__ T Box<T>::computeCircumscribedRadius() const
{
    return (norm(m_extent));
}
// -------------------------------------------------------------------------------------------------
// Returns the bounding box to box
template <typename T>
__HOSTDEVICE__ Vector3<T> Box<T>::computeBoundingBox() const
{
    return (Vector3<T>(m_extent[X], m_extent[Y], m_extent[Z]));
}
// -------------------------------------------------------------------------------------------------
// Returns the bounding cylinder to Box
template <typename T>
__HOSTDEVICE__ Vector3<T> Box<T>::computeBoundingCylinder() const
{
    // [radius, halfHeight, axisIndex(0=X,1=Y,2=Z)]
    const T ex = m_extent[X], ey = m_extent[Y], ez = m_extent[Z];
    if(ez <= ex && ez <= ey)
        return Vector3<T>(sqrt(ex * ex + ey * ey), ez, T(2));
    else if(ey <= ex && ey <= ez)
        return Vector3<T>(sqrt(ex * ex + ez * ez), ey, T(1));
    else
        return Vector3<T>(sqrt(ey * ey + ez * ez), ex, T(0));
}
// -------------------------------------------------------------------------------------------------
// Box support function, returns the support point P, i.e. the point on the
// surface of the box that satisfies max(P.v)
template <typename T>
__HOSTDEVICE__ Vector3<T> Box<T>::support(const Vector3<T>& v) const
{
    return (Vector3<T>(v[X] < T(0) ? -m_extent[X] : m_extent[X],
                       v[Y] < T(0) ? -m_extent[Y] : m_extent[Y],
                       v[Z] < T(0) ? -m_extent[Z] : m_extent[Z]));
}

#endif
