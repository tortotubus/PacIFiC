#ifndef GRAINS_GEOMETRY_SPHERE_IMPL_HH
#define GRAINS_GEOMETRY_SPHERE_IMPL_HH

#include "Sphere.hh"
#include "VectorMath.hh"

#define visuNodeNbPerQar 8

// -------------------------------------------------------------------------------------------------
// Constructor with radius
template <typename T>
__HOSTDEVICE__ Sphere<T>::Sphere(T r)
    : m_radius(r)
{
}
// -------------------------------------------------------------------------------------------------
// Destructor
template <typename T>
__HOSTDEVICE__ Sphere<T>::~Sphere()
{
}
// -------------------------------------------------------------------------------------------------
// Returns the convex type
template <typename T>
__HOSTDEVICE__ ConvexType Sphere<T>::getConvexType() const
{
    return (SPHERE);
}
// -------------------------------------------------------------------------------------------------
// Returns the radius
template <typename T>
__HOSTDEVICE__ T Sphere<T>::getRadius() const
{
    return (m_radius);
}
// -------------------------------------------------------------------------------------------------
// Sets the radius
template <typename T>
__HOSTDEVICE__ void Sphere<T>::setRadius(T r)
{
    m_radius = r;
}
// -------------------------------------------------------------------------------------------------
// Returns a clone of the sphere
template <typename T>
__HOSTDEVICE__ Convex<T>* Sphere<T>::clone() const
{
    return (new Sphere<T>(m_radius));
}
// -------------------------------------------------------------------------------------------------
// Returns the volume of the Sphere
template <typename T>
__HOSTDEVICE__ T Sphere<T>::computeVolume() const
{
    return (T(4) * PI<T> * m_radius * m_radius * m_radius / T(3));
}
// -------------------------------------------------------------------------------------------------
// Computes the diagonal inertia tensor
template <typename T>
__HOSTDEVICE__ void Sphere<T>::computeInertia(T (&inertia)[3]) const
{
    T r5 = m_radius * m_radius * m_radius * m_radius * m_radius;
    // Diagonal components: Ixx = Iyy = Izz for sphere
    inertia[2] = inertia[1] = inertia[0] = T(8) * PI<T> / T(15) * r5;
}
// -------------------------------------------------------------------------------------------------
// Returns the circumscribed radius of the Sphere
template <typename T>
__HOSTDEVICE__ T Sphere<T>::computeCircumscribedRadius() const
{
    return (m_radius);
}
// -------------------------------------------------------------------------------------------------
// Returns the bounding box to Sphere
template <typename T>
__HOSTDEVICE__ Vector3<T> Sphere<T>::computeBoundingBox() const
{
    return (Vector3<T>(m_radius, m_radius, m_radius));
}
// -------------------------------------------------------------------------------------------------
// Returns the bounding cylinder to Sphere
template <typename T>
__HOSTDEVICE__ Vector3<T> Sphere<T>::computeBoundingCylinder() const
{
    // [radius, halfHeight, axisIndex=Z(2)] -- any axis is equivalent for a sphere
    return Vector3<T>(m_radius, m_radius, T(2));
}
// -------------------------------------------------------------------------------------------------
// Sphere support function, returns the support point P, i.e. the point on the
// surface of the Sphere that satisfies max(P.v)
template <typename T>
__HOSTDEVICE__ Vector3<T> Sphere<T>::support(const Vector3<T>& v) const
{
    return (m_radius * v.normalized());
}

#endif
