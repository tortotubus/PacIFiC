#ifndef GRAINS_GEOMETRY_TRANSFORM3_IMPL_HH
#define GRAINS_GEOMETRY_TRANSFORM3_IMPL_HH

#include "Transform3.hh"
#include "GrainsUtils.hh"
#include "MatrixMath.hh"

// -------------------------------------------------------------------------------------------------
// Default constructor. Origin is def and matrix is identity
template <typename T>
__HOSTDEVICE__ Transform3<T>::Transform3(T def)
{
    m_basis.setValue(T(1), T(0), T(0), T(0), T(1), T(0), T(0), T(0), T(1));
    m_origin.setValue(def, def, def);
}
// -------------------------------------------------------------------------------------------------
// Constructor with origin coordinates as inputs and matrix is identity
template <typename T>
__HOSTDEVICE__ Transform3<T>::Transform3(T x, T y, T z)
{
    m_basis.setValue(T(1), T(0), T(0), T(0), T(1), T(0), T(0), T(0), T(1));
    m_origin.setValue(x, y, z);
}
// -------------------------------------------------------------------------------------------------
// Constructor with a 1D array of 12 values as inputs containing the
// rotation matrix coefficients following by the origin coordinates
template <typename T>
__HOSTDEVICE__ Transform3<T>::Transform3(T const* buffer)
{
    setValue(buffer);
}
// -------------------------------------------------------------------------------------------------
// Constructor with a quaternion and position
template <typename T>
__HOSTDEVICE__ Transform3<T>::Transform3(const Quaternion<T>& q, const Vector3<T>& p)
{
    m_basis  = q.toMatrix();
    m_origin = p;
}
// -------------------------------------------------------------------------------------------------
// Constructor with a two transformations such that 'this = b2w o inv(a2w) =
// b2a'
template <typename T>
__HOSTDEVICE__ Transform3<T>::Transform3(const Transform3<T>& a2w, const Transform3<T>& b2w)
    : Transform3<T>(b2w)
{
    this->relativeToTransform(a2w);
}
// -------------------------------------------------------------------------------------------------
// Copy constructor
template <typename T>
__HOSTDEVICE__ Transform3<T>::Transform3(const Transform3<T>& t)
{
    m_basis  = t.m_basis;
    m_origin = t.m_origin;
}
// -------------------------------------------------------------------------------------------------
// Destructor
template <typename T>
__HOSTDEVICE__ Transform3<T>::~Transform3()
{
}
// -------------------------------------------------------------------------------------------------
// Gets the rotation of the transformation as a quaternion
template <typename T>
__HOSTDEVICE__ Quaternion<T> Transform3<T>::getRotation() const
{
    return (Quaternion<T>(m_basis));
}
// -------------------------------------------------------------------------------------------------
// Gets the orientation of the transformation
template <typename T>
__HOSTDEVICE__ Matrix3<T> Transform3<T>::getBasis() const
{
    return (m_basis);
}
// -------------------------------------------------------------------------------------------------
// Gets the origin of the transformation
template <typename T>
__HOSTDEVICE__ Vector3<T> Transform3<T>::getOrigin() const
{
    return (m_origin);
}
// -------------------------------------------------------------------------------------------------
// Sets the transformation with an 1D array of 12 values as inputs
template <typename T>
__HOSTDEVICE__ void Transform3<T>::setValue(T const* buffer)
{
    m_basis.setValue(buffer);
    m_origin.setValue(buffer + 9);
}
// -------------------------------------------------------------------------------------------------
// Sets the matrix part of the transformation
template <typename T>
__HOSTDEVICE__ void Transform3<T>::setBasis(const Matrix3<T>& m)
{
    m_basis = m;
}
// -------------------------------------------------------------------------------------------------
// Sets the matrix part of the transformation with specified rotations around
// each principal axis
template <typename T>
__HOSTDEVICE__ void Transform3<T>::setBasis(T aX, T aY, T aZ)
{
    m_basis = Matrix3<T>(aX, aY, aZ);
}
// -------------------------------------------------------------------------------------------------
// Sets the origin of the transformation
template <typename T>
__HOSTDEVICE__ void Transform3<T>::setOrigin(const Vector3<T>& v)
{
    m_origin = v;
}
// -------------------------------------------------------------------------------------------------
// Sets the transformation to the identity
template <typename T>
__HOSTDEVICE__ void Transform3<T>::setIDentity()
{
    m_basis.setValue(T(1), T(0), T(0), T(0), T(1), T(0), T(0), T(0), T(1));
    m_origin.setValue(T(0), T(0), T(0));
}
// -------------------------------------------------------------------------------------------------
// Set the transformation to the inverse of another transformation t
template <typename T>
__HOSTDEVICE__ void Transform3<T>::setToInverseTransform(const Transform3<T>& t, bool isRotation)
{
    if(isRotation)
        m_basis = transpose(t.m_basis);
    else
        m_basis = inverse(t.m_basis);
    m_origin.setValue((-m_basis * t.m_origin).getBuffer());
}
// -------------------------------------------------------------------------------------------------
// Composition of affine transformations: this = t2 o t1 (t1 first
// followed by t2)
template <typename T>
__HOSTDEVICE__ void Transform3<T>::setToTransformsComposition(const Transform3<T>& t1,
                                                              const Transform3<T>& t2)
{
    m_basis  = t2.m_basis * t1.m_basis;
    m_origin = t2(t1.m_origin);
}
// -------------------------------------------------------------------------------------------------
// Composition with a scaling transformation: this = this o scaling
template <typename T>
__HOSTDEVICE__ void Transform3<T>::composeWithScaling(const Vector3<T>& v)
{
    scale(m_basis, v);
}
// -------------------------------------------------------------------------------------------------
// Composition on the left by a rotation described by a transform:
// this = rot o this (this first followed by rot)
// This composition leaves the origin unchanged but does not check that rot
// is indeed a rotation
template <typename T>
__HOSTDEVICE__ void Transform3<T>::composeLeftByRotation(const Transform3<T>& t)
{
    m_basis = t.m_basis * m_basis;
}
// -------------------------------------------------------------------------------------------------
// Composition on the left by a rotation described by a quaternion:
// this = rot( quaternion ) o this ( this first followed by rot( quaternion ) )
template <typename T>
__HOSTDEVICE__ void Transform3<T>::composeLeftByRotation(const Quaternion<T>& q)
{
    const T* b  = m_basis.getBuffer();
    T        qx = q[X], qy = q[Y], qz = q[Z], qw = q[W];
    T        px, py, pz, pw;

    // We compute below the matrix product M_rot(q).basis where
    // M_rot(q) is the rotation matrix corresponding to the quaternion q
    for(int i = 0; i < 3; ++i)
    {
        px            = qy * b[6 + i] - qz * b[3 + i] + qw * b[i];
        py            = qz * b[i] - qx * b[6 + i] + qw * b[3 + i];
        pz            = qx * b[3 + i] - qy * b[i] + qw * b[6 + i];
        pw            = -qx * b[i] - qy * b[3 + i] - qz * b[6 + i];
        m_basis[X][i] = qy * pz - qz * py - pw * qx + qw * px;
        m_basis[Y][i] = qz * px - qx * pz - pw * qy + qw * py;
        m_basis[Z][i] = qx * py - qy * px - pw * qz + qw * pz;
    }

    // // We compute below the matrix product M_rot(q).basis where
    // // M_rot(q) is the rotation matrix corresponding to the quaternion q
    // for(int i = 0; i < 3; ++i)
    // {
    //     px             = qy * b[6 + i] - qz * b[3 + i] + qw * b[i];
    //     py             = qz * b[i] - qx * b[6 + i] + qw * b[3 + i];
    //     pz             = qx * b[3 + i] - qy * b[i] + qw * b[6 + i];
    //     pw             = -qx * b[i] - qy * b[3 + i] - qz * b[6 + i];
    //     m_basis[i]     = qy * pz - qz * py - pw * qx + qw * px;
    //     m_basis[3 + i] = qz * px - qx * pz - pw * qy + qw * py;
    //     m_basis[6 + i] = qx * py - qy * px - pw * qz + qw * pz;
    // }

    // const Matrix3<T>& mm = m_basis;
    // T qx = q[X], qy = q[Y], qz = q[Z], qw = q[W];
    // T px, py, pz, pw;

    // // We compute below the matrix product M_rot(q).basis where
    // // M_rot(q) is the rotation matrix corresponding to the quaternion q
    // for ( int i = 0; i < 3; ++i )
    // {
    //     px = qy * mm[Z][i] - qz * mm[Y][i] + qw * mm[X][i];
    //     py = qz * mm[X][i] - qx * mm[Z][i] + qw * mm[Y][i];
    //     pz = qx * mm[Y][i] - qy * mm[X][i] + qw * mm[Z][i];
    //     pw = - qx * mm[X][i] - qy * mm[Y][i] - qz * mm[Z][i];
    //     m_basis[X][i] = qy * pz - qz * py - pw * qx + qw * px;
    //     m_basis[Y][i] = qz * px - qx * pz - pw * qy + qw * py;
    //     m_basis[Z][i] = qx * py - qy * px - pw * qz + qw * pz;
    // }
}
// -------------------------------------------------------------------------------------------------
// Composition on the left by a translation:
// this = trans(vector) o this (this first followed by trans(vector))
template <typename T>
__HOSTDEVICE__ void Transform3<T>::composeLeftByTranslation(const Vector3<T>& v)
{
    m_origin += v;
}
// -------------------------------------------------------------------------------------------------
// Composition on the left by another affine transformation:
// this = t o this (this first followed by t)
template <typename T>
__HOSTDEVICE__ void Transform3<T>::composeLeftByTransform(const Transform3<T>& t)
{
    m_origin = t.m_origin + t.m_basis * m_origin;
    m_basis  = t.m_basis * m_basis;
}
// -------------------------------------------------------------------------------------------------
// Composition on the right by another affine transformation:
// this = this o t (t first followed by this)
template <typename T>
__HOSTDEVICE__ void Transform3<T>::composeRightByTransform(const Transform3<T>& t)
{
    m_origin += m_basis * t.m_origin;
    m_basis *= t.m_basis;
}
// -------------------------------------------------------------------------------------------------
// Relative transformation with respect to t
template <typename T>
__HOSTDEVICE__ void Transform3<T>::relativeToTransform(const Transform3<T>& t)
{
    Matrix3<T> const inverseRotation = transpose(t.m_basis);
    m_basis                          = inverseRotation * m_basis;
    m_origin                         = inverseRotation * (m_origin - t.m_origin);
}
// -------------------------------------------------------------------------------------------------
// Updates the transformation with a displacement and a rotation
template <typename T>
__HOSTDEVICE__ void Transform3<T>::updateTransform(const Vector3<T>&    transMotion,
                                                   const Quaternion<T>& rotMotion)
{
    this->composeLeftByTranslation(transMotion);
    this->composeLeftByRotation(rotMotion);
}
// -------------------------------------------------------------------------------------------------
// Returns the result of applying the transformation to the input vector
template <typename T>
__HOSTDEVICE__ Vector3<T> Transform3<T>::operator()(const Vector3<T>& v) const
{
    return (Vector3<T>(m_basis * v + m_origin));
}
// -------------------------------------------------------------------------------------------------
// Equal operator to another Transform object
template <typename T>
__HOSTDEVICE__ Transform3<T>& Transform3<T>::operator=(const Transform3<T>& t)
{
    if(&t != this)
    {
        m_basis  = t.m_basis;
        m_origin = t.m_origin;
    }
    return (*this);
}

#endif
