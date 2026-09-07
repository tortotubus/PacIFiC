#ifndef _COLLISIONDETECTION_HH_
#define _COLLISIONDETECTION_HH_

#include <tuple>

#include "ContactInfo.hh"
#include "Drum.hh"
#include "GJK.hh"
#include "GrainsParameters.hh"
#include "GrainsUtils.hh"
#include "MatrixMath.hh"
#include "MiscMath.hh"
#include "OBB.hh"
#include "QuaternionMath.hh"
#include "Rectangle.hh"
#include "RigidBody.hh"
#include "Transform3.hh"
#include "Trapezoid.hh"
#include "Triangle.hh"

// =================================================================================================
/** @brief The header-only file for Rigid bodies collision detections.

    Functions for collision detection between two rigid bodies.

    @author A.Yazdani - 2024 - Construction */
// =================================================================================================
/** @name CollisionDetection Low-Level Methods */
//@{
template <typename T>
struct PanelGeometry
{
    T halfWidthBottom;
    T halfWidthTop;
    T halfHeight;
};

// -------------------------------------------------------------------------------------------------
/** @brief Helper that extracts commonly-used rigid-body properties for contact handling.
    @param rbA first rigid body
    @param rbB second rigid body */
template <typename T>
__HOSTDEVICE__ static INLINE auto getPropertiesForContact(const RigidBody<T>& rbA,
                                                          const RigidBody<T>& rbB)
{
    auto sa = rbA.getPropertiesSnapshot();
    auto sb = rbB.getPropertiesSnapshot();

    const Convex<T>* convexA     = sa.convex;
    const Convex<T>* convexB     = sb.convex;
    T                crustA      = sa.crustThickness;
    T                crustB      = sb.crustThickness;
    T                circRadiusA = sa.circumscribedRadius;
    T                circRadiusB = sb.circumscribedRadius;

    T massA          = sa.mass;
    T massB          = sb.mass;
    T invMassA       = (massA == T(0)) ? T(0) : T(1) / massA;
    T invMassB       = (massB == T(0)) ? T(0) : T(1) / massB;
    T invReducedMass = 1 / (invMassA + invMassB);

    T    invRadA       = (circRadiusA == T(0)) ? T(0) : T(1) / circRadiusA;
    T    invRadB       = (circRadiusB == T(0)) ? T(0) : T(1) / circRadiusB;
    T    averageRadius = 1 / (invRadA + invRadB);
    uint materialHash  = triangularHash(sa.material, sb.material);

    return std::make_tuple(convexA,
                           crustA,
                           circRadiusA,
                           convexB,
                           crustB,
                           circRadiusB,
                           invReducedMass,
                           averageRadius,
                           materialHash);
}

// -------------------------------------------------------------------------------------------------
/** @brief ShapeData-only variant of getPropertiesForContact. Reads all needed scalar data
    directly from ShapeData fields -- no RigidBody or virtual dispatch required.
    @param sdA pre-built ShapeData for body A
    @param sdB pre-built ShapeData for body B */
template <typename T>
__HOSTDEVICE__ static INLINE auto getPropertiesForContact(const ShapeData<T>& sdA,
                                                          const ShapeData<T>& sdB)
{
    const T    crustA         = sdA.crust;
    const T    crustB         = sdB.crust;
    const T    circRadiusA    = sdA.circumscribedRadius;
    const T    circRadiusB    = sdB.circumscribedRadius;
    const T    invReducedMass = T(1) / (sdA.invMass + sdB.invMass);
    const T    invRadA        = (circRadiusA == T(0)) ? T(0) : T(1) / circRadiusA;
    const T    invRadB        = (circRadiusB == T(0)) ? T(0) : T(1) / circRadiusB;
    const T    averageRadius  = T(1) / (invRadA + invRadB);
    const uint materialHash   = triangularHash(sdA.material, sdB.material);
    return std::make_tuple(crustA,
                           circRadiusA,
                           crustB,
                           circRadiusB,
                           invReducedMass,
                           averageRadius,
                           materialHash);
}

// -------------------------------------------------------------------------------------------------
// Returns whether two rigid bodies os spherical shape intersect
template <typename T>
__HOSTDEVICE__ static INLINE bool
    intersectSpheres(const RigidBody<T>& rbA, const RigidBody<T>& rbB, const Vector3<T>& v_b2a)
{
    T radiiSum = rbA.getCircumscribedRadius() + rbB.getCircumscribedRadius();
    T dist2    = norm2(v_b2a);
    return (dist2 < radiiSum * radiiSum);
}

// -------------------------------------------------------------------------------------------------
template <typename T>
__HOSTDEVICE__ static INLINE T clampSymmetric(T value, T limit)
{
    return value > limit ? limit : (value < -limit ? -limit : value);
}

// -------------------------------------------------------------------------------------------------
template <typename T>
__HOSTDEVICE__ static INLINE void computeDrumContactLocal(
    const Vector3<T>& cB, T radius, T halfHeight, Vector3<T>& ptA, Vector3<T>& normal)
{
    const T radial = sqrt(cB[X] * cB[X] + cB[Z] * cB[Z]);
    if(radial > LOWEPS<T>)
    {
        const T invRadial = T(1) / radial;
        ptA[X]            = radius * cB[X] * invRadial;
        ptA[Z]            = radius * cB[Z] * invRadial;
    }
    else
    {
        ptA[X] = radius;
        ptA[Z] = T(0);
    }
    ptA[Y] = clampSymmetric(cB[Y], halfHeight);

    Vector3<T> d    = cB - ptA;
    const T    dist = norm(d);
    if(dist > LOWEPS<T>)
    {
        normal = d / dist;
        return;
    }

    if(radial > LOWEPS<T>)
    {
        const T invRadial = T(1) / radial;
        const T sign      = radial >= radius ? T(1) : T(-1);
        normal            = Vector3<T>(sign * cB[X] * invRadial, T(0), sign * cB[Z] * invRadial);
        return;
    }

    normal = Vector3<T>(T(-1), T(0), T(0));
}

// -------------------------------------------------------------------------------------------------
__HOSTDEVICE__ static INLINE bool isPanelType(ConvexType type)
{
    return type == ConvexType::RECTANGLE || type == ConvexType::TRAPEZOID
           || type == ConvexType::TRIANGLE || type == ConvexType::DRUM;
}

// -------------------------------------------------------------------------------------------------
template <typename T>
__HOSTDEVICE__ static INLINE PanelGeometry<T> getPanelGeometry(const Convex<T>& convex)
{
    if(convex.getConvexType() == ConvexType::RECTANGLE)
    {
        const Vector3<T> bbox = convex.computeBoundingBox();
        return {bbox[X], bbox[X], bbox[Y]};
    }
    if(convex.getConvexType() == ConvexType::TRAPEZOID)
    {
        const auto& trap = static_cast<const Trapezoid<T>&>(convex);
        return {trap.getHalfWidthBottom(), trap.getHalfWidthTop(), trap.getHalfHeight()};
    }
    if(convex.getConvexType() == ConvexType::TRIANGLE)
    {
        const auto& tri = static_cast<const Triangle<T>&>(convex);
        return {tri.getHalfBase(), T(0), tri.getHalfHeight()};
    }
    return {T(0), T(0), T(0)};
}

// -------------------------------------------------------------------------------------------------
template <typename T>
__HOSTDEVICE__ static INLINE PanelGeometry<T> getPanelGeometry(const ShapeData<T>& shapeData)
{
    if(shapeData.type == ConvexType::RECTANGLE)
        return {shapeData.params[0], shapeData.params[0], shapeData.params[1]};
    if(shapeData.type == ConvexType::TRAPEZOID)
        return {shapeData.params[0], shapeData.params[1], shapeData.params[2]};
    if(shapeData.type == ConvexType::TRIANGLE)
        return {shapeData.params[0], T(0), shapeData.params[1]};
    return {T(0), T(0), T(0)};
}

// -------------------------------------------------------------------------------------------------
/** @brief Computes the closest point on a trapezoid panel surface to the particle center cB
    (all in the panel's body-local frame) and the outward normal from the surface toward the
    particle.
    Body-frame convention: trapezoid lies in XY plane (Z=0); bottom edge (halfWidthBottom) at
    y = -halfHeight, top edge (halfWidthTop) at y = +halfHeight.
    @param cB             particle center in panel-local frame
    @param halfWidthBottom half-width at y = -halfHeight
    @param halfWidthTop    half-width at y = +halfHeight
    @param halfHeight      half-height of the panel
    @param ptA        output: closest point on the panel surface
    @param normal     output: unit normal from ptA toward cB */
template <typename T>
__HOSTDEVICE__ static INLINE void computeTrapezoidContactLocal(const Vector3<T>& cB,
                                                               T                 halfWidthBottom,
                                                               T                 halfWidthTop,
                                                               T                 halfHeight,
                                                               Vector3<T>&       ptA,
                                                               Vector3<T>&       normal)
{
    // Clamp y to panel height
    const T cy = cB[Y] > halfHeight ? halfHeight : (cB[Y] < -halfHeight ? -halfHeight : cB[Y]);
    // Interpolate half-width at clamped y
    const T t  = (cy + halfHeight) / (T(2) * halfHeight);
    const T hw = halfWidthBottom + (halfWidthTop - halfWidthBottom) * t;
    // Clamp x to interpolated half-width
    const T cx = cB[X] > hw ? hw : (cB[X] < -hw ? -hw : cB[X]);

    ptA[X] = cx;
    ptA[Y] = cy;
    ptA[Z] = T(0);

    // Normal: from panel's closest point toward particle center
    Vector3<T> d    = cB - ptA;
    const T    dist = norm(d);
    if(dist > LOWEPS<T>)
    {
        normal = d / dist;
        return;
    }

    // Degenerate: particle centre is exactly on the panel surface; use face normal
    normal = cB[Z] >= T(0) ? Vector3<T>(T(0), T(0), T(1)) : Vector3<T>(T(0), T(0), T(-1));
}

// -------------------------------------------------------------------------------------------------
template <typename T>
__HOSTDEVICE__ static INLINE void computePanelContactLocal(const Vector3<T>& cB,
                                                           const Convex<T>&  convex,
                                                           Vector3<T>&       ptA,
                                                           Vector3<T>&       normal)
{
    if(convex.getConvexType() == ConvexType::DRUM)
    {
        const auto& drum = static_cast<const Drum<T>&>(convex);
        computeDrumContactLocal(cB, drum.getRadius(), drum.getHeight() / T(2), ptA, normal);
        return;
    }

    const auto geometry = getPanelGeometry(convex);
    computeTrapezoidContactLocal(cB,
                                 geometry.halfWidthBottom,
                                 geometry.halfWidthTop,
                                 geometry.halfHeight,
                                 ptA,
                                 normal);
}

// -------------------------------------------------------------------------------------------------
template <typename T>
__HOSTDEVICE__ static INLINE void computePanelContactLocal(const Vector3<T>&   cB,
                                                           const ShapeData<T>& shapeData,
                                                           Vector3<T>&         ptA,
                                                           Vector3<T>&         normal)
{
    if(shapeData.type == ConvexType::DRUM)
    {
        computeDrumContactLocal(cB, shapeData.params[0], shapeData.params[1], ptA, normal);
        return;
    }

    const auto geometry = getPanelGeometry(shapeData);
    computeTrapezoidContactLocal(cB,
                                 geometry.halfWidthBottom,
                                 geometry.halfWidthTop,
                                 geometry.halfHeight,
                                 ptA,
                                 normal);
}
//@}

// -------------------------------------------------------------------------------------------------
/** @brief Fills the contact-overlap fields of snapshot for the flat-panel case.
    Sets overlapDistance = normal · (ptB − ptA); fills contactPoint and contactVector
    only if overlapDistance < T(0) (actual contact).
    @param ptA closest point on panel A (world frame)
    @param ptB closest point on convex B (world frame)
    @param normal outward unit normal of panel A (world frame)
    @param snapshot snapshot to fill */
template <typename T>
__HOSTDEVICE__ static INLINE void
    fillPanelContactSnapshot(const Vector3<T>&                  ptA,
                             const Vector3<T>&                  ptB,
                             const Vector3<T>&                  normal,
                             typename ContactInfo<T>::Snapshot& snapshot)
{
    Vector3<T> n = normal;
    round(n);
    snapshot.overlapDistance = n * (ptB - ptA);
    snapshot.contactPoint    = T(0.5) * (ptA + ptB);
    snapshot.contactVector   = n;
}

// -------------------------------------------------------------------------------------------------
/** @brief Fills the contact-overlap fields of snapshot for the GJK case.
    Subtracts crustA + crustB from snapshot.overlapDistance (which must hold the raw
    GJK separation distance at call time), then fills contactPoint and contactVector.
    @param ptA world-frame witness point on body A
    @param ptB world-frame witness point on body B
    @param crustA crust thickness of body A
    @param crustB crust thickness of body B
    @param snapshot snapshot to update */
template <typename T>
__HOSTDEVICE__ static INLINE void
    fillGJKContactSnapshot(const Vector3<T>&                  ptA,
                           const Vector3<T>&                  ptB,
                           const T                            crustA,
                           const T                            crustB,
                           typename ContactInfo<T>::Snapshot& snapshot)
{
    snapshot.overlapDistance -= crustA + crustB;
    snapshot.contactPoint = T(0.5) * (ptA + ptB);
    Vector3<T> d          = (ptB - ptA).normalized();
    round(d);
    snapshot.contactVector = d;
}

/* ============================================================================================== */
/* High-Level Methods                                                                             */
/* ============================================================================================== */
/** @name CollisionDetection High-Level Methods */
//@{
/** @brief Returns whether 2 rigid bodies intersect - relative transformation.
    @param rbA first rigid body
    @param rbB second rigid body
    @param b2a geometric transformation describing convex B in the A's reference frame */
template <typename T>
__HOSTDEVICE__ inline bool
    intersectRigidBodies(const RigidBody<T>& rbA, const RigidBody<T>& rbB, const Transform3<T>& b2a)
{
    const Convex<T>& convexA = *(rbA.getConvex());
    const Convex<T>& convexB = *(rbB.getConvex());
    return (intersectGJK(convexA, convexB, b2a));
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns whether 2 rigid bodies intersect.
    @param rbA first rigid body
    @param rbB second rigid body
    @param a2w geometric transformation describing convex A in the world reference frame
    @param b2w geometric transformation describing convex B in the world reference frame */
template <typename T>
__HOSTDEVICE__ inline bool intersectRigidBodies(const RigidBody<T>&  rbA,
                                                const RigidBody<T>&  rbB,
                                                const Transform3<T>& a2w,
                                                const Transform3<T>& b2w)
{
    const Convex<T>& convexA = *(rbA.getConvex());
    const Convex<T>& convexB = *(rbB.getConvex());
    return (intersectGJK(convexA, convexB, a2w, b2w));
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns whether 2 rigid bodies intersect - relative transformation.
    @param rbA first rigid body
    @param rbB second rigid body
    @param v_b2a position describing convex B in the A's reference frame
    @param q_b2a rotation describing convex B in the A's reference frame */
template <typename T>
__HOSTDEVICE__ inline bool intersectRigidBodies(const RigidBody<T>&  rbA,
                                                const RigidBody<T>&  rbB,
                                                const Vector3<T>&    v_b2a,
                                                const Quaternion<T>& q_b2a)
{
    const Convex<T>& convexA = *(rbA.getConvex());
    const Convex<T>& convexB = *(rbB.getConvex());
    return (intersectGJK(convexA, convexB, v_b2a, q_b2a));
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns whether 2 rigid bodies intersect.
    @param rbA first rigid body
    @param rbB second rigid body
    @param v_a2w position describing convex A in the world reference frame
    @param v_b2w position describing convex B in the world reference frame
    @param q_a2w rotation describing convex A in the world reference frame
    @param q_b2w rotation describing convex B in the world reference frame */
template <typename T>
__HOSTDEVICE__ inline bool intersectRigidBodies(const RigidBody<T>&  rbA,
                                                const RigidBody<T>&  rbB,
                                                const Vector3<T>&    v_a2w,
                                                const Vector3<T>&    v_b2w,
                                                const Quaternion<T>& q_a2w,
                                                const Quaternion<T>& q_b2w)
{
    const Convex<T>& convexA = *(rbA.getConvex());
    const Convex<T>& convexB = *(rbB.getConvex());
    return (intersectGJK(convexA, convexB, v_a2w, v_b2w, q_a2w, q_b2w));
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns the contact information (if any) for 2 rigid bodies - relative transformation.
    @param rbA first rigid body
    @param rbB second rigid body
    @param b2a geometric transformation describing convex B in the A's reference frame
    @param contactInfo output contact information */
template <typename T,
          GJKType            GJKVARIANT,
          bool               GJKACC,
          BoundingVolumeType BVType = BoundingVolumeType::OFF>
__HOSTDEVICE__ inline void closestPointsRigidBodies(const RigidBody<T>&  rbA,
                                                    const RigidBody<T>&  rbB,
                                                    const Transform3<T>& b2a,
                                                    ContactInfo<T>&      contactInfo)
{
    // BV pre-filter: sphere reject then OBB SAT before entering GJK
    if constexpr(BVType == BoundingVolumeType::OBB)
    {
        const T radiiSum = rbA.getCircumscribedRadius() + rbB.getCircumscribedRadius();
        if(norm2(b2a.getOrigin()) >= radiiSum * radiiSum)
        {
            contactInfo.setOverlapDistance(T(1));
            return;
        }
        if(!intersectOrientedBoundingBox(rbA.getConvex()->computeBoundingBox(),
                                         rbB.getConvex()->computeBoundingBox(),
                                         b2a))
        {
            contactInfo.setOverlapDistance(T(1));
            return;
        }
    }
    /* ---------------------------------------------------------------------------------------------
    Comments on the contactInfo. It applies to all variants of this function:
    1. If actual overlap distance (GJK dist - crustA - crustB < 0), there is contact otherwise no
    contact. Although we can enforce an early exit, but other threads are most likely still
    running, so we continue to have consistent code path.

    2. ptA and ptB are in their respective local coordinate systems and represent points on the
    actual rigid bodies, not the shrunken versions. Contact point definition as the mid point
    between ptA and ptB.

    3. If contact, overlap is negative and overlap_vector is from B to A If no contact, overlap is
    positive and we do not care about the direction of overlap_vector. Assuming A and B are the
    centers of the 2 convex bodies overlap_vector = overlap * Vector3(A to B)
    --------------------------------------------------------------------------------------------- */
    // Extract properties for contact and populate snapshot
    auto [convexA_ptr,
          crustA,
          circRadiusA,
          convexB_ptr,
          crustB,
          circRadiusB,
          averageMass,
          averageRadius,
          contactHash]
        = getPropertiesForContact(rbA, rbB);
    typename ContactInfo<T>::Snapshot snapshot;
    snapshot.averageMass     = averageMass;
    snapshot.averageRadius   = averageRadius;
    snapshot.contactHash     = contactHash;
    snapshot.overlapDistance = std::numeric_limits<T>::max();

    // Get convexes and their types
    const Convex<T>& convexA = *(convexA_ptr);
    const Convex<T>& convexB = *(convexB_ptr);
    const ConvexType typeA   = convexA.getConvexType();
    const ConvexType typeB   = convexB.getConvexType();

    // Sphere-Sphere Case
    if(typeA == ConvexType::SPHERE && typeB == ConvexType::SPHERE)
    {
        T rA                     = rbA.getCircumscribedRadius();
        T rB                     = rbB.getCircumscribedRadius();
        snapshot.contactVector   = b2a.getOrigin();
        snapshot.overlapDistance = norm(snapshot.contactVector) - rA - rB;
        snapshot.contactPoint    = (rA + T(.5) * snapshot.overlapDistance) * snapshot.contactVector;
        snapshot.contactVector.normalize();
    }
    // Flat-panel Case: rectangle / trapezoid / triangle via the same local closest-point rule.
    else if(isPanelType(typeA))
    {
        Vector3<T> ptA, ptB, normal;
        computePanelContactLocal(b2a.getOrigin(), convexA, ptA, normal);
        ptB = (b2a)(convexB.support((-normal) * b2a.getBasis()));
        fillPanelContactSnapshot(ptA, ptB, normal, snapshot);
    }
    else if(isPanelType(typeB))
    {
        GAbort("General wall collision detection is not implemented yet.");
    }
    // General Case
    else if constexpr(GJKVARIANT != GJKType::SPHERE)
    {
        Vector3<T> ptA, ptB;
        uint       nbIterGJK     = 0;
        snapshot.overlapDistance = computeClosestPoints_GJK<T, GJKVARIANT, GJKACC>(convexA,
                                                                                   convexB,
                                                                                   b2a,
                                                                                   crustA,
                                                                                   crustB,
                                                                                   ptA,
                                                                                   ptB,
                                                                                   nbIterGJK);
        if(snapshot.overlapDistance < HIGHEPS<T>)
        {
            // GJK degenerate (deep overlap): retry with inflated crust (factor 10, min 50%
            // of shape size) to recover a contact direction, then impose nominal overlap.
            const Vector3<T> bboxA       = convexA.computeBoundingBox();
            const Vector3<T> bboxB       = convexB.computeBoundingBox();
            const T          inscribedRA = min(min(bboxA[X], bboxA[Y]), bboxA[Z]);
            const T          inscribedRB = min(min(bboxB[X], bboxB[Y]), bboxB[Z]);
            const T          retryCA     = min(T(10) * crustA, inscribedRA);
            const T          retryCB     = min(T(10) * crustB, inscribedRB);
            Vector3<T>       rptA, rptB;
            uint             rIter     = 0;
            T                retryDist = computeClosestPoints_GJK<T, GJKVARIANT, GJKACC>(convexA,
                                                                          convexB,
                                                                          b2a,
                                                                          retryCA,
                                                                          retryCB,
                                                                          rptA,
                                                                          rptB,
                                                                          rIter);
            if(retryDist >= HIGHEPS<T>)
            {
                rptB                   = (b2a)(rptB);
                snapshot.contactPoint  = T(0.5) * (rptA + rptB);
                snapshot.contactVector = (rptB - rptA).normalized();
            }
            else
            {
                snapshot.contactVector = b2a.getOrigin().normalized();
                snapshot.contactPoint  = T(0.5) * b2a.getOrigin();
            }
            snapshot.overlapDistance = -(crustA + crustB);
        }
        else
        {
            // ptA = (a2a)(ptA);
            ptB = (b2a)(ptB);
            fillGJKContactSnapshot(ptA, ptB, crustA, crustB, snapshot);
        }
    }

    // Set contact information
    contactInfo.setSnapshot(snapshot);
    return;
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns the contact information (if any) for 2 rigid bodies.
    @param rbA first rigid body
    @param rbB second rigid body
    @param a2w geometric transformation describing convex A in the world reference frame
    @param b2w geometric transformation describing convex B in the world reference frame
    @param contactInfo output contact information */
template <typename T,
          GJKType            GJKVARIANT,
          bool               GJKACC,
          BoundingVolumeType BVType = BoundingVolumeType::OFF>
__HOSTDEVICE__ inline void closestPointsRigidBodies(const RigidBody<T>&  rbA,
                                                    const RigidBody<T>&  rbB,
                                                    const Transform3<T>& a2w,
                                                    const Transform3<T>& b2w,
                                                    ContactInfo<T>&      contactInfo)
{
    // BV pre-filter: sphere reject then OBB or OBC SAT before entering GJK
    if constexpr(BVType == BoundingVolumeType::OBB)
    {
        const T radiiSum = rbA.getCircumscribedRadius() + rbB.getCircumscribedRadius();
        if(norm2(b2w.getOrigin() - a2w.getOrigin()) >= radiiSum * radiiSum)
        {
            contactInfo.setOverlapDistance(T(1));
            return;
        }
        if(!intersectOrientedBoundingBox(rbA.getConvex()->computeBoundingBox(),
                                         rbB.getConvex()->computeBoundingBox(),
                                         a2w,
                                         b2w))
        {
            contactInfo.setOverlapDistance(T(1));
            return;
        }
    }
    if constexpr(BVType == BoundingVolumeType::OBC)
    {
        const T radiiSum = rbA.getCircumscribedRadius() + rbB.getCircumscribedRadius();
        if(norm2(b2w.getOrigin() - a2w.getOrigin()) >= radiiSum * radiiSum)
        {
            contactInfo.setOverlapDistance(T(1));
            return;
        }
        auto axisFromIndex = [](T idx) -> Vector3<T> {
            if(idx == T(0))
                return Vector3<T>(T(1), T(0), T(0));
            if(idx == T(1))
                return Vector3<T>(T(0), T(1), T(0));
            return Vector3<T>(T(0), T(0), T(1));
        };
        const Vector3<T> bcA = rbA.getConvex()->computeBoundingCylinder();
        const Vector3<T> bcB = rbB.getConvex()->computeBoundingCylinder();
        if(!intersectOrientedBoundingCylinder(bcA[X],
                                              bcA[Y],
                                              axisFromIndex(bcA[Z]),
                                              bcB[X],
                                              bcB[Y],
                                              axisFromIndex(bcB[Z]),
                                              a2w,
                                              b2w))
        {
            contactInfo.setOverlapDistance(T(1));
            return;
        }
    }
    // Extract properties for contact and populate snapshot
    auto [convexA_ptr,
          crustA,
          circRadiusA,
          convexB_ptr,
          crustB,
          circRadiusB,
          averageMass,
          averageRadius,
          contactHash]
        = getPropertiesForContact(rbA, rbB);
    typename ContactInfo<T>::Snapshot snapshot;
    snapshot.averageMass     = averageMass;
    snapshot.averageRadius   = averageRadius;
    snapshot.contactHash     = contactHash;
    snapshot.overlapDistance = std::numeric_limits<T>::max();

    // Get convexes and their types
    const Convex<T>& convexA = *(convexA_ptr);
    const Convex<T>& convexB = *(convexB_ptr);
    const ConvexType typeA   = convexA.getConvexType();
    const ConvexType typeB   = convexB.getConvexType();

    Vector3<T> ptA, ptB;

    // Sphere-Sphere Case
    if(typeA == ConvexType::SPHERE && typeB == ConvexType::SPHERE)
    {
        T rA                     = rbA.getCircumscribedRadius();
        T rB                     = rbB.getCircumscribedRadius();
        ptA                      = a2w.getOrigin();
        ptB                      = b2w.getOrigin() - ptA;
        snapshot.overlapDistance = norm(ptB) - rA - rB;
        {
            Vector3<T> d_hat       = ptB.normalized();
            snapshot.contactPoint  = ptA + (rA + T(.5) * snapshot.overlapDistance) * d_hat;
            snapshot.contactVector = d_hat;
        }
    }
    // Flat-panel Case: rectangle / trapezoid / triangle via the same local closest-point rule.
    else if(isPanelType(typeA))
    {
        const Vector3<T>& vA       = a2w.getOrigin();
        const Matrix3<T>& mA       = a2w.getBasis();
        const Vector3<T>  cB_local = inverse(mA) * (b2w.getOrigin() - vA);
        Vector3<T>        ptA_local, normal_local;
        computePanelContactLocal(cB_local, convexA, ptA_local, normal_local);
        const Vector3<T> normal = mA * normal_local;
        ptA                     = vA + mA * ptA_local;
        ptB                     = (b2w)(convexB.support((-normal) * b2w.getBasis()));
        fillPanelContactSnapshot(ptA, ptB, normal, snapshot);
    }
    else if(isPanelType(typeB))
    {
        GAbort("General wall collision detection is not implemented yet.");
    }
    // General Case
    else if constexpr(GJKVARIANT != GJKType::SPHERE)
    {
        uint nbIterGJK           = 0;
        snapshot.overlapDistance = computeClosestPoints_GJK<T, GJKVARIANT, GJKACC>(convexA,
                                                                                   convexB,
                                                                                   a2w,
                                                                                   b2w,
                                                                                   crustA,
                                                                                   crustB,
                                                                                   ptA,
                                                                                   ptB,
                                                                                   nbIterGJK);
        if(snapshot.overlapDistance < HIGHEPS<T>)
        {
            const Vector3<T> bboxA       = convexA.computeBoundingBox();
            const Vector3<T> bboxB       = convexB.computeBoundingBox();
            const T          inscribedRA = min(min(bboxA[X], bboxA[Y]), bboxA[Z]);
            const T          inscribedRB = min(min(bboxB[X], bboxB[Y]), bboxB[Z]);
            const T          retryCA     = min(T(10) * crustA, inscribedRA);
            const T          retryCB     = min(T(10) * crustB, inscribedRB);
            Vector3<T>       rptA, rptB;
            uint             rIter     = 0;
            T                retryDist = computeClosestPoints_GJK<T, GJKVARIANT, GJKACC>(convexA,
                                                                          convexB,
                                                                          a2w,
                                                                          b2w,
                                                                          retryCA,
                                                                          retryCB,
                                                                          rptA,
                                                                          rptB,
                                                                          rIter);
            if(retryDist >= HIGHEPS<T>)
            {
                rptA                   = (a2w)(rptA);
                rptB                   = (b2w)(rptB);
                snapshot.contactPoint  = T(0.5) * (rptA + rptB);
                snapshot.contactVector = (rptB - rptA).normalized();
            }
            else
            {
                snapshot.contactVector = (b2w.getOrigin() - a2w.getOrigin()).normalized();
                snapshot.contactPoint  = T(0.5) * (a2w.getOrigin() + b2w.getOrigin());
            }
            snapshot.overlapDistance = -(crustA + crustB);
        }
        else
        {
            ptA = (a2w)(ptA);
            ptB = (b2w)(ptB);
            fillGJKContactSnapshot(ptA, ptB, crustA, crustB, snapshot);
        }
    }

    // Set contact information
    contactInfo.setSnapshot(snapshot);
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns the contact information (if any) for 2 rigid bodies - relative transformation.
    @param rbA first rigid body
    @param rbB second rigid body
    @param v_b2a position describing convex B in the A's reference frame
    @param q_b2a rotation describing convex B in the A's reference frame
    @param contactInfo output contact information */
template <typename T,
          GJKType            GJKVARIANT,
          bool               GJKACC,
          BoundingVolumeType BVType = BoundingVolumeType::OFF>
__HOSTDEVICE__ inline void closestPointsRigidBodies(const RigidBody<T>&  rbA,
                                                    const RigidBody<T>&  rbB,
                                                    const Vector3<T>&    v_b2a,
                                                    const Quaternion<T>& q_b2a,
                                                    ContactInfo<T>&      contactInfo)
{
    // BV pre-filter: sphere reject then OBB SAT before entering GJK
    if constexpr(BVType == BoundingVolumeType::OBB)
    {
        const T radiiSum = rbA.getCircumscribedRadius() + rbB.getCircumscribedRadius();
        if(norm2(v_b2a) >= radiiSum * radiiSum)
        {
            contactInfo.setOverlapDistance(T(1));
            return;
        }
        if(!intersectOrientedBoundingBox(rbA.getConvex()->computeBoundingBox(),
                                         rbB.getConvex()->computeBoundingBox(),
                                         v_b2a,
                                         q_b2a))
        {
            contactInfo.setOverlapDistance(T(1));
            return;
        }
    }
    // Extract properties for contact and populate snapshot
    auto [convexA_ptr,
          crustA,
          circRadiusA,
          convexB_ptr,
          crustB,
          circRadiusB,
          averageMass,
          averageRadius,
          contactHash]
        = getPropertiesForContact(rbA, rbB);
    typename ContactInfo<T>::Snapshot snapshot;
    snapshot.averageMass     = averageMass;
    snapshot.averageRadius   = averageRadius;
    snapshot.contactHash     = contactHash;
    snapshot.overlapDistance = std::numeric_limits<T>::max();

    // Get convexes and their types
    const Convex<T>& convexA = *(convexA_ptr);
    const Convex<T>& convexB = *(convexB_ptr);
    const ConvexType typeA   = convexA.getConvexType();
    const ConvexType typeB   = convexB.getConvexType();

    // Sphere-Sphere Case
    if(typeA == ConvexType::SPHERE && typeB == ConvexType::SPHERE)
    {
        T rA                     = rbA.getCircumscribedRadius();
        T rB                     = rbB.getCircumscribedRadius();
        snapshot.contactVector   = v_b2a;
        snapshot.overlapDistance = norm(snapshot.contactVector) - rA - rB;
        snapshot.contactPoint    = (rA + T(.5) * snapshot.overlapDistance) * snapshot.contactVector;
        snapshot.contactVector.normalize();
    }
    // Flat-panel Case: rectangle / trapezoid / triangle via the same local closest-point rule.
    else if(isPanelType(typeA))
    {
        Vector3<T> ptA, ptB, normal;
        computePanelContactLocal(v_b2a, convexA, ptA, normal);
        ptB = convexB.support(q_b2a << (-normal));
        transform(q_b2a, v_b2a, ptB);
        fillPanelContactSnapshot(ptA, ptB, normal, snapshot);
    }
    else if(isPanelType(typeB))
    {
        GAbort("General wall collision detection is not implemented yet.");
    }
    // General Case
    else if constexpr(GJKVARIANT != GJKType::SPHERE)
    {
        Vector3<T> ptA, ptB;
        uint       nbIterGJK     = 0;
        snapshot.overlapDistance = computeClosestPoints_GJK<T, GJKVARIANT, GJKACC>(convexA,
                                                                                   convexB,
                                                                                   v_b2a,
                                                                                   q_b2a,
                                                                                   crustA,
                                                                                   crustB,
                                                                                   ptA,
                                                                                   ptB,
                                                                                   nbIterGJK);
        if(snapshot.overlapDistance < HIGHEPS<T>)
        {
            const Vector3<T> bboxA       = convexA.computeBoundingBox();
            const Vector3<T> bboxB       = convexB.computeBoundingBox();
            const T          inscribedRA = min(min(bboxA[X], bboxA[Y]), bboxA[Z]);
            const T          inscribedRB = min(min(bboxB[X], bboxB[Y]), bboxB[Z]);
            const T          retryCA     = min(T(10) * crustA, inscribedRA);
            const T          retryCB     = min(T(10) * crustB, inscribedRB);
            Vector3<T>       rptA, rptB;
            uint             rIter     = 0;
            T                retryDist = computeClosestPoints_GJK<T, GJKVARIANT, GJKACC>(convexA,
                                                                          convexB,
                                                                          v_b2a,
                                                                          q_b2a,
                                                                          retryCA,
                                                                          retryCB,
                                                                          rptA,
                                                                          rptB,
                                                                          rIter);
            if(retryDist >= HIGHEPS<T>)
            {
                transform(q_b2a, v_b2a, rptB);
                snapshot.contactPoint  = T(0.5) * (rptA + rptB);
                snapshot.contactVector = (rptB - rptA).normalized();
            }
            else
            {
                snapshot.contactVector = v_b2a.normalized();
                snapshot.contactPoint  = T(0.5) * v_b2a;
            }
            snapshot.overlapDistance = -(crustA + crustB);
        }
        else
        {
            // transform(q_a2a, v_a2a, ptA);
            transform(q_b2a, v_b2a, ptB);
            fillGJKContactSnapshot(ptA, ptB, crustA, crustB, snapshot);
        }
    }

    // Set contact information
    contactInfo.setSnapshot(snapshot);
    return;
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns the contact information (if any) for 2 rigid bodies.
    @param rbA first rigid body
    @param rbB second rigid body
    @param v_a2w position describing convex A in the world reference frame
    @param v_b2w position describing convex B in the world reference frame
    @param q_a2w rotation describing convex A in the world reference frame
    @param q_b2w rotation describing convex B in the world reference frame
    @param contactInfo output contact information */
template <typename T,
          GJKType            GJKVARIANT,
          bool               GJKACC,
          BoundingVolumeType BVType = BoundingVolumeType::OFF>
__HOSTDEVICE__ inline void closestPointsRigidBodies(const RigidBody<T>&  rbA,
                                                    const RigidBody<T>&  rbB,
                                                    const Vector3<T>&    v_a2w,
                                                    const Vector3<T>&    v_b2w,
                                                    const Quaternion<T>& q_a2w,
                                                    const Quaternion<T>& q_b2w,
                                                    ContactInfo<T>&      contactInfo)
{
    // BV pre-filter: sphere reject then OBB or OBC SAT before entering GJK
    if constexpr(BVType == BoundingVolumeType::OBB)
    {
        const T radiiSum = rbA.getCircumscribedRadius() + rbB.getCircumscribedRadius();
        if(norm2(v_b2w - v_a2w) >= radiiSum * radiiSum)
        {
            contactInfo.setOverlapDistance(T(1));
            return;
        }
        if(!intersectOrientedBoundingBox(rbA.getConvex()->computeBoundingBox(),
                                         rbB.getConvex()->computeBoundingBox(),
                                         v_a2w,
                                         v_b2w,
                                         q_a2w,
                                         q_b2w))
        {
            contactInfo.setOverlapDistance(T(1));
            return;
        }
    }
    if constexpr(BVType == BoundingVolumeType::OBC)
    {
        const T radiiSum = rbA.getCircumscribedRadius() + rbB.getCircumscribedRadius();
        if(norm2(v_b2w - v_a2w) >= radiiSum * radiiSum)
        {
            contactInfo.setOverlapDistance(T(1));
            return;
        }
        auto axisFromIndex = [](T idx) -> Vector3<T> {
            if(idx == T(0))
                return Vector3<T>(T(1), T(0), T(0));
            if(idx == T(1))
                return Vector3<T>(T(0), T(1), T(0));
            return Vector3<T>(T(0), T(0), T(1));
        };
        const Vector3<T> bcA = rbA.getConvex()->computeBoundingCylinder();
        const Vector3<T> bcB = rbB.getConvex()->computeBoundingCylinder();
        if(!intersectOrientedBoundingCylinder(bcA[X],
                                              bcA[Y],
                                              axisFromIndex(bcA[Z]),
                                              bcB[X],
                                              bcB[Y],
                                              axisFromIndex(bcB[Z]),
                                              v_a2w,
                                              v_b2w,
                                              q_a2w,
                                              q_b2w))
        {
            contactInfo.setOverlapDistance(T(1));
            return;
        }
    }
    // Extract properties for contact and populate snapshot
    auto [convexA_ptr,
          crustA,
          circRadiusA,
          convexB_ptr,
          crustB,
          circRadiusB,
          averageMass,
          averageRadius,
          contactHash]
        = getPropertiesForContact(rbA, rbB);
    typename ContactInfo<T>::Snapshot snapshot;
    snapshot.averageMass     = averageMass;
    snapshot.averageRadius   = averageRadius;
    snapshot.contactHash     = contactHash;
    snapshot.overlapDistance = std::numeric_limits<T>::max();

    // Get convexes and their types
    const Convex<T>& convexA = *(convexA_ptr);
    const Convex<T>& convexB = *(convexB_ptr);
    const ConvexType typeA   = convexA.getConvexType();
    const ConvexType typeB   = convexB.getConvexType();

    Vector3<T> ptA, ptB;

    // Sphere-Sphere Case
    if(typeA == ConvexType::SPHERE && typeB == ConvexType::SPHERE)
    {
        T rA                     = rbA.getCircumscribedRadius();
        T rB                     = rbB.getCircumscribedRadius();
        ptA                      = v_a2w;
        ptB                      = v_b2w - ptA;
        snapshot.overlapDistance = norm(ptB) - rA - rB;
        {
            Vector3<T> d_hat       = ptB.normalized();
            snapshot.contactPoint  = ptA + (rA + T(.5) * snapshot.overlapDistance) * d_hat;
            snapshot.contactVector = d_hat;
        }
    }
    // Flat-panel Case: rectangle / trapezoid / triangle via the same local closest-point rule.
    else if(isPanelType(typeA))
    {
        const Vector3<T> cB_local = q_a2w << (v_b2w - v_a2w);
        Vector3<T>       ptA_local, normal_local;
        computePanelContactLocal(cB_local, convexA, ptA_local, normal_local);
        const Vector3<T> normal = q_a2w >> normal_local;
        ptA                     = v_a2w + (q_a2w >> ptA_local);
        ptB                     = (q_b2w >> convexB.support(q_b2w << (-normal))) + v_b2w;
        fillPanelContactSnapshot(ptA, ptB, normal, snapshot);
    }
    else if(isPanelType(typeB))
    {
        GAbort("General wall collision detection is not implemented yet.");
    }
    // General Case
    else if constexpr(GJKVARIANT != GJKType::SPHERE)
    {
        uint nbIterGJK           = 0;
        snapshot.overlapDistance = computeClosestPoints_GJK<T, GJKVARIANT, GJKACC>(convexA,
                                                                                   convexB,
                                                                                   v_a2w,
                                                                                   v_b2w,
                                                                                   q_a2w,
                                                                                   q_b2w,
                                                                                   crustA,
                                                                                   crustB,
                                                                                   ptA,
                                                                                   ptB,
                                                                                   nbIterGJK);
        if(snapshot.overlapDistance < HIGHEPS<T>)
        {
            const Vector3<T> bboxA       = convexA.computeBoundingBox();
            const Vector3<T> bboxB       = convexB.computeBoundingBox();
            const T          inscribedRA = min(min(bboxA[X], bboxA[Y]), bboxA[Z]);
            const T          inscribedRB = min(min(bboxB[X], bboxB[Y]), bboxB[Z]);
            const T          retryCA     = min(T(10) * crustA, inscribedRA);
            const T          retryCB     = min(T(10) * crustB, inscribedRB);
            Vector3<T>       rptA, rptB;
            uint             rIter     = 0;
            T                retryDist = computeClosestPoints_GJK<T, GJKVARIANT, GJKACC>(convexA,
                                                                          convexB,
                                                                          v_a2w,
                                                                          v_b2w,
                                                                          q_a2w,
                                                                          q_b2w,
                                                                          retryCA,
                                                                          retryCB,
                                                                          rptA,
                                                                          rptB,
                                                                          rIter);
            if(retryDist >= HIGHEPS<T>)
            {
                transform(q_a2w, v_a2w, rptA);
                transform(q_b2w, v_b2w, rptB);
                snapshot.contactPoint  = T(0.5) * (rptA + rptB);
                snapshot.contactVector = (rptB - rptA).normalized();
            }
            else
            {
                snapshot.contactVector = (v_b2w - v_a2w).normalized();
                snapshot.contactPoint  = T(0.5) * (v_a2w + v_b2w);
            }
            snapshot.overlapDistance = -(crustA + crustB);
        }
        else
        {
            transform(q_a2w, v_a2w, ptA);
            transform(q_b2w, v_b2w, ptB);
            fillGJKContactSnapshot(ptA, ptB, crustA, crustB, snapshot);
        }
    }

    // Set contact information
    contactInfo.setSnapshot(snapshot);
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns the distance between 2 rigid bodies.
    @param rbA first rigid body
    @param rbB second rigid body
    @param a2w geometric transformation describing convex A in the world reference frame
    @param b2w geometric transformation describing convex B in the world reference frame
    @param method method identifier (currently unused) */
template <typename T, GJKType GJKVARIANT, bool GJKACC>
__HOSTDEVICE__ inline T distanceRigidBodies(const RigidBody<T>&  rbA,
                                            const RigidBody<T>&  rbB,
                                            const Transform3<T>& a2w,
                                            const Transform3<T>& b2w,
                                            const uint           method)
{
    Convex<T> const* convexA = rbA.getConvex();
    Convex<T> const* convexB = rbB.getConvex();

    Vector3<T> ptA, ptB;
    uint       nbIterGJK = 0;
    T          distance  = 0;
    distance             = computeClosestPoints_GJK<T, GJKVARIANT, GJKACC>(*convexA,
                                                               *convexB,
                                                               a2w,
                                                               b2w,
                                                               rbA.getCrustThickness(),
                                                               rbB.getCrustThickness(),
                                                               ptA,
                                                               ptB,
                                                               nbIterGJK);
    return (distance);
}

// -------------------------------------------------------------------------------------------------
/** @brief RigidBody-free version: all contact properties read from ShapeData directly.
    Used by the prebuilt GPU path (BVType always OFF; BV filter ran separately).
    Rectangle-Particle: replaces virtual convexB.support with device_support_raw and
    isInside with direct params[] bounds check.
    @param sdA pre-built ShapeData for body A (shapeId-indexed)
    @param sdB pre-built ShapeData for body B (shapeId-indexed)
    @param v_b2a position of B in A-local frame
    @param q_b2a rotation of B in A-local frame
    @param contactInfo output contact information */
template <typename T, GJKType GJKVARIANT, bool GJKACC>
__HOSTDEVICE__ inline void closestPointsRigidBodies(const ShapeData<T>&  sdA,
                                                    const ShapeData<T>&  sdB,
                                                    const Vector3<T>&    v_b2a,
                                                    const Quaternion<T>& q_b2a,
                                                    ContactInfo<T>&      contactInfo)
{
    auto [crustA, circRadiusA, crustB, circRadiusB, averageMass, averageRadius, contactHash]
        = getPropertiesForContact(sdA, sdB);

    typename ContactInfo<T>::Snapshot snapshot;
    snapshot.averageMass     = averageMass;
    snapshot.averageRadius   = averageRadius;
    snapshot.contactHash     = contactHash;
    snapshot.overlapDistance = std::numeric_limits<T>::max();

    const ConvexType typeA = sdA.type;
    const ConvexType typeB = sdB.type;

    // Sphere-Sphere Case
    if(typeA == ConvexType::SPHERE && typeB == ConvexType::SPHERE)
    {
        snapshot.contactVector   = v_b2a;
        snapshot.overlapDistance = norm(snapshot.contactVector) - circRadiusA - circRadiusB;
        snapshot.contactPoint
            = (circRadiusA + T(.5) * snapshot.overlapDistance) * snapshot.contactVector;
        snapshot.contactVector.normalize();
    }
    // Flat-panel Case: rectangle / trapezoid / triangle via the same local closest-point rule.
    else if(isPanelType(typeA))
    {
        Vector3<T> ptA, ptB, normal;
        computePanelContactLocal(v_b2a, sdA, ptA, normal);
        ptB = device_support_raw(sdB, q_b2a << (-normal));
        transform(q_b2a, v_b2a, ptB);
        fillPanelContactSnapshot(ptA, ptB, normal, snapshot);
    }
    else if(isPanelType(typeB))
    {
        GAbort("General wall collision detection is not implemented yet.");
    }
    // General Case: vtable-free GJK via ShapeData
    else if constexpr(GJKVARIANT != GJKType::SPHERE)
    {
        Vector3<T> ptA, ptB;
        uint       nbIterGJK     = 0;
        snapshot.overlapDistance = computeClosestPoints_GJK<T, GJKVARIANT, GJKACC>(sdA,
                                                                                   sdB,
                                                                                   v_b2a,
                                                                                   q_b2a,
                                                                                   crustA,
                                                                                   crustB,
                                                                                   ptA,
                                                                                   ptB,
                                                                                   nbIterGJK);
        if(snapshot.overlapDistance < HIGHEPS<T>)
        {
            const T    retryCA = min(T(10) * crustA, sdA.inscribedRadius);
            const T    retryCB = min(T(10) * crustB, sdB.inscribedRadius);
            Vector3<T> rptA, rptB;
            uint       rIter     = 0;
            T          retryDist = computeClosestPoints_GJK<T, GJKVARIANT, GJKACC>(sdA,
                                                                          sdB,
                                                                          v_b2a,
                                                                          q_b2a,
                                                                          retryCA,
                                                                          retryCB,
                                                                          rptA,
                                                                          rptB,
                                                                          rIter);
            if(retryDist >= HIGHEPS<T>)
            {
                transform(q_b2a, v_b2a, rptB);
                snapshot.contactPoint  = T(0.5) * (rptA + rptB);
                snapshot.contactVector = (rptB - rptA).normalized();
            }
            else
            {
                snapshot.contactVector = v_b2a.normalized();
                snapshot.contactPoint  = T(0.5) * v_b2a;
            }
            snapshot.overlapDistance = -(crustA + crustB);
        }
        else
        {
            transform(q_b2a, v_b2a, ptB);
            fillGJKContactSnapshot(ptA, ptB, crustA, crustB, snapshot);
        }
    }

    contactInfo.setSnapshot(snapshot);
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns the contact information (if any) for 2 rigid bodies using pre-built ShapeData
    for vtable-free GJK support evaluation -- relative transformation (vec/quat).
    The General Case replaces virtual dispatch with device_support() over ShapeData.
    Sphere-Sphere and Rectangle-Particle paths fall through to the standard virtual variants.
    @param rbA first rigid body
    @param rbB second rigid body
    @param sdA pre-built ShapeData for rbA (slot-indexed)
    @param sdB pre-built ShapeData for rbB (slot-indexed)
    @param v_b2a position describing convex B in the A's reference frame
    @param q_b2a rotation describing convex B in the A's reference frame
    @param contactInfo output contact information */
template <typename T,
          GJKType            GJKVARIANT,
          bool               GJKACC,
          BoundingVolumeType BVType = BoundingVolumeType::OFF>
__HOSTDEVICE__ inline void closestPointsRigidBodies(const RigidBody<T>&  rbA,
                                                    const RigidBody<T>&  rbB,
                                                    const ShapeData<T>&  sdA,
                                                    const ShapeData<T>&  sdB,
                                                    const Vector3<T>&    v_b2a,
                                                    const Quaternion<T>& q_b2a,
                                                    ContactInfo<T>&      contactInfo)
{
    // BV pre-filter: sphere reject then OBB SAT before entering GJK
    if constexpr(BVType == BoundingVolumeType::OBB)
    {
        const T radiiSum = rbA.getCircumscribedRadius() + rbB.getCircumscribedRadius();
        if(norm2(v_b2a) >= radiiSum * radiiSum)
        {
            contactInfo.setOverlapDistance(T(1));
            return;
        }
        if(!intersectOrientedBoundingBox(rbA.getConvex()->computeBoundingBox(),
                                         rbB.getConvex()->computeBoundingBox(),
                                         v_b2a,
                                         q_b2a))
        {
            contactInfo.setOverlapDistance(T(1));
            return;
        }
    }
    // Extract properties for contact and populate snapshot
    auto [convexA_ptr,
          crustA,
          circRadiusA,
          convexB_ptr,
          crustB,
          circRadiusB,
          averageMass,
          averageRadius,
          contactHash]
        = getPropertiesForContact(rbA, rbB);
    typename ContactInfo<T>::Snapshot snapshot;
    snapshot.averageMass     = averageMass;
    snapshot.averageRadius   = averageRadius;
    snapshot.contactHash     = contactHash;
    snapshot.overlapDistance = std::numeric_limits<T>::max();

    // Get convexes and their types
    const Convex<T>& convexA = *(convexA_ptr);
    const Convex<T>& convexB = *(convexB_ptr);
    const ConvexType typeA   = convexA.getConvexType();
    const ConvexType typeB   = convexB.getConvexType();

    // Sphere-Sphere Case
    if(typeA == ConvexType::SPHERE && typeB == ConvexType::SPHERE)
    {
        T rA                     = rbA.getCircumscribedRadius();
        T rB                     = rbB.getCircumscribedRadius();
        snapshot.contactVector   = v_b2a;
        snapshot.overlapDistance = norm(snapshot.contactVector) - rA - rB;
        snapshot.contactPoint    = (rA + T(.5) * snapshot.overlapDistance) * snapshot.contactVector;
        snapshot.contactVector.normalize();
    }
    // Flat-panel Case: rectangle / trapezoid / triangle via the same local closest-point rule.
    else if(isPanelType(typeA))
    {
        Vector3<T> ptA, ptB, normal;
        computePanelContactLocal(v_b2a, convexA, ptA, normal);
        ptB = convexB.support(q_b2a << (-normal));
        transform(q_b2a, v_b2a, ptB);
        fillPanelContactSnapshot(ptA, ptB, normal, snapshot);
    }
    else if(isPanelType(typeB))
    {
        GAbort("General wall collision detection is not implemented yet.");
    }
    // General Case: vtable-free via ShapeData
    else if constexpr(GJKVARIANT != GJKType::SPHERE)
    {
        Vector3<T> ptA, ptB;
        uint       nbIterGJK     = 0;
        snapshot.overlapDistance = computeClosestPoints_GJK<T, GJKVARIANT, GJKACC>(sdA,
                                                                                   sdB,
                                                                                   v_b2a,
                                                                                   q_b2a,
                                                                                   crustA,
                                                                                   crustB,
                                                                                   ptA,
                                                                                   ptB,
                                                                                   nbIterGJK);
        if(snapshot.overlapDistance < HIGHEPS<T>)
        {
            const T    retryCA = min(T(10) * crustA, sdA.inscribedRadius);
            const T    retryCB = min(T(10) * crustB, sdB.inscribedRadius);
            Vector3<T> rptA, rptB;
            uint       rIter     = 0;
            T          retryDist = computeClosestPoints_GJK<T, GJKVARIANT, GJKACC>(sdA,
                                                                          sdB,
                                                                          v_b2a,
                                                                          q_b2a,
                                                                          retryCA,
                                                                          retryCB,
                                                                          rptA,
                                                                          rptB,
                                                                          rIter);
            if(retryDist >= HIGHEPS<T>)
            {
                transform(q_b2a, v_b2a, rptB);
                snapshot.contactPoint  = T(0.5) * (rptA + rptB);
                snapshot.contactVector = (rptB - rptA).normalized();
            }
            else
            {
                snapshot.contactVector = v_b2a.normalized();
                snapshot.contactPoint  = T(0.5) * v_b2a;
            }
            snapshot.overlapDistance = -(crustA + crustB);
        }
        else
        {
            // transform(q_a2a, v_a2a, ptA);
            transform(q_b2a, v_b2a, ptB);
            fillGJKContactSnapshot(ptA, ptB, crustA, crustB, snapshot);
        }
    }

    // Set contact information
    contactInfo.setSnapshot(snapshot);
    return;
}

// -------------------------------------------------------------------------------------------------
/** @brief RigidBody-free world-frame version: all contact properties from ShapeData.
    Used by the prebuilt GPU path (BVType always OFF; BV filter ran separately).
    @param sdA pre-built ShapeData for body A (shapeId-indexed)
    @param sdB pre-built ShapeData for body B (shapeId-indexed)
    @param v_a2w world-frame position of A
    @param v_b2w world-frame position of B
    @param q_a2w world-frame orientation of A
    @param q_b2w world-frame orientation of B
    @param contactInfo output contact information */
template <typename T, GJKType GJKVARIANT, bool GJKACC>
__HOSTDEVICE__ inline void closestPointsRigidBodies(const ShapeData<T>&  sdA,
                                                    const ShapeData<T>&  sdB,
                                                    const Vector3<T>&    v_a2w,
                                                    const Vector3<T>&    v_b2w,
                                                    const Quaternion<T>& q_a2w,
                                                    const Quaternion<T>& q_b2w,
                                                    ContactInfo<T>&      contactInfo)
{
    auto [crustA, circRadiusA, crustB, circRadiusB, averageMass, averageRadius, contactHash]
        = getPropertiesForContact(sdA, sdB);

    typename ContactInfo<T>::Snapshot snapshot;
    snapshot.averageMass     = averageMass;
    snapshot.averageRadius   = averageRadius;
    snapshot.contactHash     = contactHash;
    snapshot.overlapDistance = std::numeric_limits<T>::max();

    const ConvexType typeA = sdA.type;
    const ConvexType typeB = sdB.type;
    Vector3<T>       ptA, ptB;

    // Sphere-Sphere Case
    if(typeA == ConvexType::SPHERE && typeB == ConvexType::SPHERE)
    {
        ptA                      = v_a2w;
        ptB                      = v_b2w - ptA;
        snapshot.overlapDistance = norm(ptB) - circRadiusA - circRadiusB;
        {
            Vector3<T> d_hat       = ptB.normalized();
            snapshot.contactPoint  = ptA + (circRadiusA + T(.5) * snapshot.overlapDistance) * d_hat;
            snapshot.contactVector = d_hat;
        }
    }
    // Flat-panel Case: rectangle / trapezoid / triangle via the same local closest-point rule.
    else if(isPanelType(typeA))
    {
        const Vector3<T> cB_local = q_a2w << (v_b2w - v_a2w);
        Vector3<T>       ptA_local, normal_local;
        computePanelContactLocal(cB_local, sdA, ptA_local, normal_local);
        const Vector3<T> normal = q_a2w >> normal_local;
        ptA                     = v_a2w + (q_a2w >> ptA_local);
        ptB                     = (q_b2w >> device_support_raw(sdB, q_b2w << (-normal))) + v_b2w;
        fillPanelContactSnapshot(ptA, ptB, normal, snapshot);
    }
    else if(isPanelType(typeB))
    {
        GAbort("General wall collision detection is not implemented yet.");
    }
    // General Case: vtable-free GJK via ShapeData
    else if constexpr(GJKVARIANT != GJKType::SPHERE)
    {
        uint nbIterGJK           = 0;
        snapshot.overlapDistance = computeClosestPoints_GJK<T, GJKVARIANT, GJKACC>(sdA,
                                                                                   sdB,
                                                                                   v_a2w,
                                                                                   v_b2w,
                                                                                   q_a2w,
                                                                                   q_b2w,
                                                                                   crustA,
                                                                                   crustB,
                                                                                   ptA,
                                                                                   ptB,
                                                                                   nbIterGJK);
        if(snapshot.overlapDistance < HIGHEPS<T>)
        {
            const T    retryCA = min(T(10) * crustA, sdA.inscribedRadius);
            const T    retryCB = min(T(10) * crustB, sdB.inscribedRadius);
            Vector3<T> rptA, rptB;
            uint       rIter     = 0;
            T          retryDist = computeClosestPoints_GJK<T, GJKVARIANT, GJKACC>(sdA,
                                                                          sdB,
                                                                          v_a2w,
                                                                          v_b2w,
                                                                          q_a2w,
                                                                          q_b2w,
                                                                          retryCA,
                                                                          retryCB,
                                                                          rptA,
                                                                          rptB,
                                                                          rIter);
            if(retryDist >= HIGHEPS<T>)
            {
                transform(q_a2w, v_a2w, rptA);
                transform(q_b2w, v_b2w, rptB);
                snapshot.contactPoint  = T(0.5) * (rptA + rptB);
                snapshot.contactVector = (rptB - rptA).normalized();
            }
            else
            {
                snapshot.contactVector = (v_b2w - v_a2w).normalized();
                snapshot.contactPoint  = T(0.5) * (v_a2w + v_b2w);
            }
            snapshot.overlapDistance = -(crustA + crustB);
        }
        else
        {
            transform(q_a2w, v_a2w, ptA);
            transform(q_b2w, v_b2w, ptB);
            fillGJKContactSnapshot(ptA, ptB, crustA, crustB, snapshot);
        }
    }

    contactInfo.setSnapshot(snapshot);
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns the contact information (if any) for 2 rigid bodies using pre-built ShapeData
    for vtable-free GJK support evaluation -- world-frame (vec/quat).
    The General Case replaces virtual dispatch with device_support() over ShapeData.
    Sphere-Sphere and Rectangle-Particle paths fall through to the standard virtual variants.
    @param rbA first rigid body
    @param rbB second rigid body
    @param sdA pre-built ShapeData for rbA (slot-indexed)
    @param sdB pre-built ShapeData for rbB (slot-indexed)
    @param v_a2w position describing convex A in the world reference frame
    @param v_b2w position describing convex B in the world reference frame
    @param q_a2w rotation describing convex A in the world reference frame
    @param q_b2w rotation describing convex B in the world reference frame
    @param contactInfo output contact information */
template <typename T,
          GJKType            GJKVARIANT,
          bool               GJKACC,
          BoundingVolumeType BVType = BoundingVolumeType::OFF>
__HOSTDEVICE__ inline void closestPointsRigidBodies(const RigidBody<T>&  rbA,
                                                    const RigidBody<T>&  rbB,
                                                    const ShapeData<T>&  sdA,
                                                    const ShapeData<T>&  sdB,
                                                    const Vector3<T>&    v_a2w,
                                                    const Vector3<T>&    v_b2w,
                                                    const Quaternion<T>& q_a2w,
                                                    const Quaternion<T>& q_b2w,
                                                    ContactInfo<T>&      contactInfo)
{
    // BV pre-filter: sphere reject then OBB or OBC SAT before entering GJK
    if constexpr(BVType == BoundingVolumeType::OBB)
    {
        const T radiiSum = rbA.getCircumscribedRadius() + rbB.getCircumscribedRadius();
        if(norm2(v_b2w - v_a2w) >= radiiSum * radiiSum)
        {
            contactInfo.setOverlapDistance(T(1));
            return;
        }
        if(!intersectOrientedBoundingBox(rbA.getConvex()->computeBoundingBox(),
                                         rbB.getConvex()->computeBoundingBox(),
                                         v_a2w,
                                         v_b2w,
                                         q_a2w,
                                         q_b2w))
        {
            contactInfo.setOverlapDistance(T(1));
            return;
        }
    }
    if constexpr(BVType == BoundingVolumeType::OBC)
    {
        const T radiiSum = rbA.getCircumscribedRadius() + rbB.getCircumscribedRadius();
        if(norm2(v_b2w - v_a2w) >= radiiSum * radiiSum)
        {
            contactInfo.setOverlapDistance(T(1));
            return;
        }
        auto axisFromIndex = [](T idx) -> Vector3<T> {
            if(idx == T(0))
                return Vector3<T>(T(1), T(0), T(0));
            if(idx == T(1))
                return Vector3<T>(T(0), T(1), T(0));
            return Vector3<T>(T(0), T(0), T(1));
        };
        const Vector3<T> bcA = rbA.getConvex()->computeBoundingCylinder();
        const Vector3<T> bcB = rbB.getConvex()->computeBoundingCylinder();
        if(!intersectOrientedBoundingCylinder(bcA[X],
                                              bcA[Y],
                                              axisFromIndex(bcA[Z]),
                                              bcB[X],
                                              bcB[Y],
                                              axisFromIndex(bcB[Z]),
                                              v_a2w,
                                              v_b2w,
                                              q_a2w,
                                              q_b2w))
        {
            contactInfo.setOverlapDistance(T(1));
            return;
        }
    }
    // Extract properties for contact and populate snapshot
    auto [convexA_ptr,
          crustA,
          circRadiusA,
          convexB_ptr,
          crustB,
          circRadiusB,
          averageMass,
          averageRadius,
          contactHash]
        = getPropertiesForContact(rbA, rbB);
    typename ContactInfo<T>::Snapshot snapshot;
    snapshot.averageMass     = averageMass;
    snapshot.averageRadius   = averageRadius;
    snapshot.contactHash     = contactHash;
    snapshot.overlapDistance = std::numeric_limits<T>::max();

    // Get convexes and their types
    const Convex<T>& convexA = *(convexA_ptr);
    const Convex<T>& convexB = *(convexB_ptr);
    const ConvexType typeA   = convexA.getConvexType();
    const ConvexType typeB   = convexB.getConvexType();

    Vector3<T> ptA, ptB;

    // Sphere-Sphere Case
    if(typeA == ConvexType::SPHERE && typeB == ConvexType::SPHERE)
    {
        T rA                     = rbA.getCircumscribedRadius();
        T rB                     = rbB.getCircumscribedRadius();
        ptA                      = v_a2w;
        ptB                      = v_b2w - ptA;
        snapshot.overlapDistance = norm(ptB) - rA - rB;
        {
            Vector3<T> d_hat       = ptB.normalized();
            snapshot.contactPoint  = ptA + (rA + T(.5) * snapshot.overlapDistance) * d_hat;
            snapshot.contactVector = d_hat;
        }
    }
    // Flat-panel Case: rectangle / trapezoid / triangle via the same local closest-point rule.
    else if(isPanelType(typeA))
    {
        const Vector3<T> cB_local = q_a2w << (v_b2w - v_a2w);
        Vector3<T>       ptA_local, normal_local;
        computePanelContactLocal(cB_local, convexA, ptA_local, normal_local);
        const Vector3<T> normal = q_a2w >> normal_local;
        ptA                     = v_a2w + (q_a2w >> ptA_local);
        ptB                     = (q_b2w >> convexB.support(q_b2w << (-normal))) + v_b2w;
        fillPanelContactSnapshot(ptA, ptB, normal, snapshot);
    }
    else if(isPanelType(typeB))
    {
        GAbort("General wall collision detection is not implemented yet.");
    }
    // General Case: vtable-free via ShapeData
    else if constexpr(GJKVARIANT != GJKType::SPHERE)
    {
        uint nbIterGJK           = 0;
        snapshot.overlapDistance = computeClosestPoints_GJK<T, GJKVARIANT, GJKACC>(sdA,
                                                                                   sdB,
                                                                                   v_a2w,
                                                                                   v_b2w,
                                                                                   q_a2w,
                                                                                   q_b2w,
                                                                                   crustA,
                                                                                   crustB,
                                                                                   ptA,
                                                                                   ptB,
                                                                                   nbIterGJK);
        if(snapshot.overlapDistance < HIGHEPS<T>)
        {
            const T    retryCA = min(T(10) * crustA, sdA.inscribedRadius);
            const T    retryCB = min(T(10) * crustB, sdB.inscribedRadius);
            Vector3<T> rptA, rptB;
            uint       rIter     = 0;
            T          retryDist = computeClosestPoints_GJK<T, GJKVARIANT, GJKACC>(sdA,
                                                                          sdB,
                                                                          v_a2w,
                                                                          v_b2w,
                                                                          q_a2w,
                                                                          q_b2w,
                                                                          retryCA,
                                                                          retryCB,
                                                                          rptA,
                                                                          rptB,
                                                                          rIter);
            if(retryDist >= HIGHEPS<T>)
            {
                transform(q_a2w, v_a2w, rptA);
                transform(q_b2w, v_b2w, rptB);
                snapshot.contactPoint  = T(0.5) * (rptA + rptB);
                snapshot.contactVector = (rptB - rptA).normalized();
            }
            else
            {
                snapshot.contactVector = (v_b2w - v_a2w).normalized();
                snapshot.contactPoint  = T(0.5) * (v_a2w + v_b2w);
            }
            snapshot.overlapDistance = -(crustA + crustB);
        }
        else
        {
            transform(q_a2w, v_a2w, ptA);
            transform(q_b2w, v_b2w, ptB);
            fillGJKContactSnapshot(ptA, ptB, crustA, crustB, snapshot);
        }
    }

    // Set contact information
    contactInfo.setSnapshot(snapshot);
}
//@}

#endif
