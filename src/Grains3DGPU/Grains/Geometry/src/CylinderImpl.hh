#ifndef GRAINS_GEOMETRY_CYLINDER_IMPL_HH
#define GRAINS_GEOMETRY_CYLINDER_IMPL_HH

#include "Cylinder.hh"
#include "VectorMath.hh"

// multiple of 4
#define visuNodeNbOnPer 32

// -------------------------------------------------------------------------------------------------
// Constructor with radius and height as input parameters
template <typename T>
__HOSTDEVICE__ Cylinder<T>::Cylinder(T r, T h)
    : m_radius(r)
    , m_halfHeight(h / T(2))
{
}
// -------------------------------------------------------------------------------------------------
// Destructor
template <typename T>
__HOSTDEVICE__ Cylinder<T>::~Cylinder()
{
}
// -------------------------------------------------------------------------------------------------
// Returns the convex type
template <typename T>
__HOSTDEVICE__ ConvexType Cylinder<T>::getConvexType() const
{
    return (CYLINDER);
}
// -------------------------------------------------------------------------------------------------
// Returns the radius
template <typename T>
__HOSTDEVICE__ T Cylinder<T>::getRadius() const
{
    return (m_radius);
}
// -------------------------------------------------------------------------------------------------
// Returns the height
template <typename T>
__HOSTDEVICE__ T Cylinder<T>::getHeight() const
{
    return (T(2) * m_halfHeight);
}
// -------------------------------------------------------------------------------------------------
// Returns a clone of the cylinder
template <typename T>
__HOSTDEVICE__ Convex<T>* Cylinder<T>::clone() const
{
    return (new Cylinder<T>(m_radius, T(2) * m_halfHeight));
}
// -------------------------------------------------------------------------------------------------
// Returns the volume of the cylinder
template <typename T>
__HOSTDEVICE__ T Cylinder<T>::computeVolume() const
{
    return (TWO_PI<T> * m_halfHeight * m_radius * m_radius);
}
// -------------------------------------------------------------------------------------------------
// Computes the diagonal inertia tensor
template <typename T>
__HOSTDEVICE__ void Cylinder<T>::computeInertia(T (&inertia)[3]) const
{
    T r2 = m_radius * m_radius;
    T c  = T(.5) * PI<T> * m_halfHeight * r2;
    // Diagonal components: Ixx = Izz (perpendicular to axis), Iyy (along axis)
    inertia[0] = inertia[2] = c * (T(4) * m_halfHeight * m_halfHeight / T(3) + r2);
    inertia[1]              = T(2) * c * r2;
}
// -------------------------------------------------------------------------------------------------
// Returns the circumscribed radius of the cylinder
template <typename T>
__HOSTDEVICE__ T Cylinder<T>::computeCircumscribedRadius() const
{

    return (sqrt(m_radius * m_radius + m_halfHeight * m_halfHeight));
}
// -------------------------------------------------------------------------------------------------
// Returns the bounding box to cylinder
template <typename T>
__HOSTDEVICE__ Vector3<T> Cylinder<T>::computeBoundingBox() const
{
    return (Vector3<T>(m_radius, m_halfHeight, m_radius));
}
// -------------------------------------------------------------------------------------------------
// Returns the bounding cylinder to Cylinder (exact fit)
template <typename T>
__HOSTDEVICE__ Vector3<T> Cylinder<T>::computeBoundingCylinder() const
{
    // [radius, halfHeight, axisIndex=Y(1)]
    return Vector3<T>(m_radius, m_halfHeight, T(1));
}
// -------------------------------------------------------------------------------------------------
// Cylinder support function, returns the support point P, i.e. the point on
// the surface of the Cylinder that satisfies max(P.v)
template <typename T>
__HOSTDEVICE__ Vector3<T> Cylinder<T>::support(const Vector3<T>& v) const
{
    T norm = sqrt(v[X] * v[X] + v[Y] * v[Y] + v[Z] * v[Z]);
    if(norm > EPS<T>)
    {
        T s = sqrt(v[X] * v[X] + v[Z] * v[Z]);
        if(s > EPS<T>)
        {
            T d = m_radius / s;
            if(fabs(v[Y]) < EPS<T>)
                return (Vector3<T>(v[X] * d, T(0), v[Z] * d));
            else
                return (Vector3<T>(v[X] * d, v[Y] < T(0) ? -m_halfHeight : m_halfHeight, v[Z] * d));
        }
        else
            return (Vector3<T>(T(0), v[Y] < T(0) ? -m_halfHeight : m_halfHeight, T(0)));
    }
    else
        return (Vector3<T>());
}

#endif
