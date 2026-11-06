#ifndef GRAINS_GEOMETRY_CONVEX_IMPL_HH
#define GRAINS_GEOMETRY_CONVEX_IMPL_HH

#include "Convex.hh"
#include "VectorMath.hh"

// -------------------------------------------------------------------------------------------------
// Default constructor
template <typename T>
__HOSTDEVICE__ Convex<T>::Convex()
{
}
// -------------------------------------------------------------------------------------------------
// Destructor
template <typename T>
__HOSTDEVICE__ Convex<T>::~Convex()
{
}
// -------------------------------------------------------------------------------------------------
// Eroded convex support function
template <typename T>
__HOSTDEVICE__ Vector3<T> Convex<T>::support(const Vector3<T>& v, T crust, T invNorm) const
{
    return support(v) - (crust * invNorm) * v;
}
// -------------------------------------------------------------------------------------------------
// Returns whether point p lies in the convex shape
// @param p point
template <typename T>
__HOSTDEVICE__ bool Convex<T>::isInside(const Vector3<T>& p) const
{
    // Default implementation (for convex shapes that are not defined)
    return true;
}

#endif
