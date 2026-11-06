#ifndef GRAINS_GEOMETRY_CONTACT_INFO_IMPL_HH
#define GRAINS_GEOMETRY_CONTACT_INFO_IMPL_HH

#include "ContactInfo.hh"
#include "GrainsUtils.hh"
#include "VectorMath.hh"

// -------------------------------------------------------------------------------------------------
// Default constructor
template <typename T>
__HOSTDEVICE__ ContactInfo<T>::ContactInfo()
    : m_contactPoint()
    , m_contactVector()
    , m_overlapDistance(T(0))
    , m_averageMass(T(0))
    , m_averageRadius(T(0))
    , m_contactHash(0)
{
}
// -------------------------------------------------------------------------------------------------
// Constructor with contact point, contact vector, and overlap distance
template <typename T>
__HOSTDEVICE__ ContactInfo<T>::ContactInfo(const Vector3<T>& pt, const Vector3<T>& vec, T overlap)
    : m_contactPoint(pt)
    , m_contactVector(vec)
    , m_overlapDistance(overlap)
    , m_averageMass(T(0))
    , m_averageRadius(T(0))
    , m_contactHash(0)
{
}
// -------------------------------------------------------------------------------------------------
// Destructor
template <typename T>
__HOSTDEVICE__ ContactInfo<T>::~ContactInfo()
{
}
// -------------------------------------------------------------------------------------------------
// Gets the contact point
template <typename T>
__HOSTDEVICE__ Vector3<T> ContactInfo<T>::getContactPoint() const
{
    return (m_contactPoint);
}
// -------------------------------------------------------------------------------------------------
// Gets the contact vector
template <typename T>
__HOSTDEVICE__ Vector3<T> ContactInfo<T>::getContactVector() const
{
    return (m_contactVector);
}
// -------------------------------------------------------------------------------------------------
// Gets the overlap distance
template <typename T>
__HOSTDEVICE__ T ContactInfo<T>::getOverlapDistance() const
{
    return (m_overlapDistance);
}
// -------------------------------------------------------------------------------------------------
// Gets the average mass
template <typename T>
__HOSTDEVICE__ T ContactInfo<T>::getAverageMass() const
{
    return m_averageMass;
}
// -------------------------------------------------------------------------------------------------
// Gets the average radius
template <typename T>
__HOSTDEVICE__ T ContactInfo<T>::getAverageRadius() const
{
    return m_averageRadius;
}
// -------------------------------------------------------------------------------------------------
// Gets the contact hash
template <typename T>
__HOSTDEVICE__ uint ContactInfo<T>::getContactHash() const
{
    return m_contactHash;
}
// -------------------------------------------------------------------------------------------------
// Gets snapshot
template <typename T>
__HOSTDEVICE__ typename ContactInfo<T>::Snapshot ContactInfo<T>::getSnapshot() const
{
    Snapshot s;
    s.contactPoint    = m_contactPoint;
    s.contactVector   = m_contactVector;
    s.overlapDistance = m_overlapDistance;
    s.averageMass     = m_averageMass;
    s.averageRadius   = m_averageRadius;
    s.contactHash     = m_contactHash;
    return s;
}
// -------------------------------------------------------------------------------------------------
// Sets the contact point
template <typename T>
__HOSTDEVICE__ void ContactInfo<T>::setContactPoint(const Vector3<T>& p)
{
    m_contactPoint = p;
}
// -------------------------------------------------------------------------------------------------
// Sets the contact vector
template <typename T>
__HOSTDEVICE__ void ContactInfo<T>::setContactVector(const Vector3<T>& v)
{
    m_contactVector = v;
}
// -------------------------------------------------------------------------------------------------
// Sets the overlap distance
template <typename T>
__HOSTDEVICE__ void ContactInfo<T>::setOverlapDistance(T d)
{
    m_overlapDistance = d;
}
// -------------------------------------------------------------------------------------------------
// Sets the average mass
template <typename T>
__HOSTDEVICE__ void ContactInfo<T>::setAverageMass(T avgMass)
{
    m_averageMass = avgMass;
}
// -------------------------------------------------------------------------------------------------
// Sets the average radius
template <typename T>
__HOSTDEVICE__ void ContactInfo<T>::setAverageRadius(T avgRadius)
{
    m_averageRadius = avgRadius;
}
// -------------------------------------------------------------------------------------------------
// Sets the contact hash
template <typename T>
__HOSTDEVICE__ void ContactInfo<T>::setContactHash(uint hash)
{
    m_contactHash = hash;
}
// -------------------------------------------------------------------------------------------------
template <typename T>
__HOSTDEVICE__ void ContactInfo<T>::setSnapshot(const ContactInfo<T>::Snapshot& s)
{
    m_contactPoint    = s.contactPoint;
    m_contactVector   = s.contactVector;
    m_overlapDistance = s.overlapDistance;
    m_averageMass     = s.averageMass;
    m_averageRadius   = s.averageRadius;
    m_contactHash     = s.contactHash;
}
// -------------------------------------------------------------------------------------------------
// Equality operator
template <typename T>
__HOSTDEVICE__ bool ContactInfo<T>::operator==(const ContactInfo<T>& other) const
{
    return (m_contactPoint == other.m_contactPoint) && (m_contactVector == other.m_contactVector)
           && (m_overlapDistance == other.m_overlapDistance)
           && (m_averageMass == other.m_averageMass) && (m_averageRadius == other.m_averageRadius)
           && (m_contactHash == other.m_contactHash);
}
// -------------------------------------------------------------------------------------------------
// Inequality operator
template <typename T>
__HOSTDEVICE__ bool ContactInfo<T>::operator!=(const ContactInfo<T>& other) const
{
    return !(*this == other);
}

#endif
