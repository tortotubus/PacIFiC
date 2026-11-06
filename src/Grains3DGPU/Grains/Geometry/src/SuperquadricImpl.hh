#ifndef GRAINS_GEOMETRY_SUPERQUADRIC_IMPL_HH
#define GRAINS_GEOMETRY_SUPERQUADRIC_IMPL_HH

#include "Superquadric.hh"
#include "MiscMath.hh"
#include "VectorMath.hh"

// multiple of 4
#define visuNodeNbOnPer 16

// -------------------------------------------------------------------------------------------------
// Constructor with half edge length as input parameters
template <typename T>
__HOSTDEVICE__ Superquadric<T>::Superquadric(T a, T b, T c, T n1, T n2)
    : m_a(a)
    , m_b(b)
    , m_c(c)
    , m_n1(n1)
    , m_n2(n2)
{
}
// -------------------------------------------------------------------------------------------------
// Destructor
template <typename T>
__HOSTDEVICE__ Superquadric<T>::~Superquadric()
{
}
// -------------------------------------------------------------------------------------------------
// Returns the convex type
template <typename T>
__HOSTDEVICE__ ConvexType Superquadric<T>::getConvexType() const
{
    return (SUPERQUADRIC);
}
// -------------------------------------------------------------------------------------------------
// Returns the extent in a Vector3 format
template <typename T>
__HOSTDEVICE__ Vector3<T> Superquadric<T>::getExtent() const
{
    return (Vector3<T>(m_a, m_b, m_c));
}
// -------------------------------------------------------------------------------------------------
// Returns the exponents (blockiness) in a Vector3 format with Z = 0
template <typename T>
__HOSTDEVICE__ Vector3<T> Superquadric<T>::getExponent() const
{
    return (Vector3<T>(m_n1, m_n2, T(0)));
}
// -------------------------------------------------------------------------------------------------
// Returns a clone of the superquadric
template <typename T>
__HOSTDEVICE__ Convex<T>* Superquadric<T>::clone() const
{
    return (new Superquadric<T>(m_a, m_b, m_c, m_n1, m_n2));
}
// -------------------------------------------------------------------------------------------------
// Returns the volume of the Superquadric
template <typename T>
__HOSTDEVICE__ T Superquadric<T>::computeVolume() const
{
    T const C    = T(2) / T(3) * m_a * m_b * m_c;
    T const eps1 = T(2) / m_n1;
    T const eps2 = T(2) / m_n2;

    /* Refer to MiscMath.hh for info about grainsBeta( x, y ). */
    return (C * eps1 * eps2 * grainsBeta(eps1, T(0.5) * eps1)
            * grainsBeta(T(0.5) * eps2, T(0.5) * eps2));
}
// -------------------------------------------------------------------------------------------------
// Computes the diagonal inertia tensor
template <typename T>
__HOSTDEVICE__ void Superquadric<T>::computeInertia(T (&inertia)[3]) const
{
    T const eps1 = T(2) / m_n1;
    T const eps2 = T(2) / m_n2;
    T const C    = T(0.4) * m_a * m_b * m_c * eps1 * eps2;

    /* Refer to MiscMath.hh for info about grainsBeta( x, y ). */
    T const prod1
        = grainsBeta(T(1.5) * eps2, T(0.5) * eps2) * grainsBeta(T(2) * eps1, T(0.5) * eps1);
    T const prod2
        = m_c * m_c * grainsBeta(T(0.5) * eps2, T(0.5) * eps2) * grainsBeta(T(1.5) * eps1, eps1);

    // Diagonal components: Ixx, Iyy, Izz
    inertia[0] = C * (m_b * m_b * prod1 + prod2);
    inertia[1] = C * (m_a * m_a * prod1 + prod2);
    inertia[2] = C * (m_a * m_a + m_b * m_b) * prod1;
}
// -------------------------------------------------------------------------------------------------
// Returns the circumscribed radius of the Superquadric
template <typename T>
__HOSTDEVICE__ T Superquadric<T>::computeCircumscribedRadius() const
{
    if((m_n1 == T(2)) && (m_n2 == T(2)))
        return (max(m_a, max(m_b, m_c)));
    else if(m_n1 == T(2))
    {
        T const alpha = pow(m_b / m_a, T(2) / (m_n2 - T(2)));
        T const xt    = T(1) / pow(T(1) + pow(alpha, m_n2), T(1) / m_n2);
        return (max(m_c, sqrt(m_a * xt * m_a * xt + alpha * m_b * xt * alpha * m_b * xt)));
    }
    else if(m_n2 == T(2))
    {
        T const m    = max(m_a, m_b);
        T const beta = pow(m_c * m_c / (m * m), T(1) / (m_n1 - T(2)));
        T const xt   = T(1) / pow(T(1) + pow(beta, m_n1), T(1) / m_n1);
        return (sqrt(m * xt * m * xt + beta * m_c * xt * beta * m_c * xt));
    }
    else
    {
        T const alpha = pow(m_b / m_a, T(2) / (m_n2 - T(2)));
        T const gamma = pow(T(1) + pow(alpha, m_n2), m_n1 / m_n2 - T(1));
        T const beta  = pow(gamma * m_c * m_c / (m_a * m_a), T(1) / (m_n1 - T(1)));
        T const xt
            = T(1) / pow(pow(T(1) + pow(alpha, m_n2), m_n1 / m_n2) + pow(beta, m_n1), T(1) / m_n1);

        return (sqrt(m_a * xt * m_a * xt + alpha * m_b * xt * alpha * m_b * xt
                     + beta * m_c * xt * beta * m_c * xt));
    }
}
// -------------------------------------------------------------------------------------------------
// Returns the bounding box to Superquadric
template <typename T>
__HOSTDEVICE__ Vector3<T> Superquadric<T>::computeBoundingBox() const
{
    return (Vector3<T>(m_a, m_b, m_c));
}
// -------------------------------------------------------------------------------------------------
// Returns the bounding cylinder to Superquadric
template <typename T>
__HOSTDEVICE__ Vector3<T> Superquadric<T>::computeBoundingCylinder() const
{
    // [radius, halfHeight, axisIndex(0=X,1=Y,2=Z)]
    if(m_c <= m_a && m_c <= m_b)
        return Vector3<T>(sqrt(m_a * m_a + m_b * m_b), m_c, T(2));
    else if(m_b <= m_a && m_b <= m_c)
        return Vector3<T>(sqrt(m_a * m_a + m_c * m_c), m_b, T(1));
    else
        return Vector3<T>(sqrt(m_b * m_b + m_c * m_c), m_a, T(0));
}
// -------------------------------------------------------------------------------------------------
// Superquadric support function, returns the support point P, i.e. the point on
// the surface of the Superquadric that satisfies max(P.v)
template <typename T>
__HOSTDEVICE__ Vector3<T> Superquadric<T>::support(const Vector3<T>& v) const
{
    T const abvx  = fabs(v[X]);
    T const abvy  = fabs(v[Y]);
    T const abvz  = fabs(v[Z]);
    T const signx = sgn(v[X]);
    T const signy = sgn(v[Y]);
    T const signz = sgn(v[Z]);

    Vector3<T> sup;

    if(abvx == T(0))
    {
        if(abvy == T(0))
            return (Vector3<T>(T(0), T(0), signz * m_c));
        else
        {
            T const alpha = pow(m_c / m_b * abvz / abvy, T(1) / (m_n1 - T(1)));
            T const yt    = T(1) / pow(T(1) + pow(alpha, m_n1), T(1) / m_n1);
            return (Vector3<T>(T(0), signy * m_b * yt, signz * alpha * m_c * yt));
        }
    }
    else
    {
        T const alpha = pow(m_b / m_a * abvy / abvx, T(1) / (m_n2 - T(1)));
        T const temp  = T(1) + pow(alpha, m_n2);
        T const gamma = pow(temp, (m_n1 - m_n2) / (m_n2 * (m_n1 - T(1))));
        T const beta  = gamma * pow(m_c / m_a * abvz / abvx, T(1) / (m_n1 - T(1)));
        T const xt    = T(1) / pow(pow(temp, m_n1 / m_n2) + pow(beta, m_n1), T(1) / m_n1);
        return (Vector3<T>(signx * m_a * xt, signy * alpha * m_b * xt, signz * beta * m_c * xt));
    }
}

#endif
