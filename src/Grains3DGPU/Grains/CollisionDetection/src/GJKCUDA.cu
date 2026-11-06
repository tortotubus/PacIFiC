#include "GJKImpl.hh"

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
#define X(T)                                                               \
    template __HOSTDEVICE__ bool intersectGJK(const Convex<T>&     a,      \
                                              const Convex<T>&     b,      \
                                              const Transform3<T>& b2a);   \
    template __HOSTDEVICE__ bool intersectGJK(const Convex<T>&     a,      \
                                              const Convex<T>&     b,      \
                                              const Transform3<T>& a2w,    \
                                              const Transform3<T>& b2w);   \
    template __HOSTDEVICE__ bool intersectGJK(const Convex<T>&     a,      \
                                              const Convex<T>&     b,      \
                                              const Vector3<T>&    v_b2a,  \
                                              const Quaternion<T>& q_b2a); \
    template __HOSTDEVICE__ bool intersectGJK(const Convex<T>&     a,      \
                                              const Convex<T>&     b,      \
                                              const Vector3<T>&    v_a2w,  \
                                              const Vector3<T>&    v_b2w,  \
                                              const Quaternion<T>& q_a2w,  \
                                              const Quaternion<T>& q_b2w);
X(float)
X(double)
#undef X

#define X(T, GJK, ACC)                                                                           \
    template __HOSTDEVICE__ T computeClosestPoints_GJK<T, GJK, ACC>(const Convex<T>&     a,      \
                                                                    const Convex<T>&     b,      \
                                                                    const Transform3<T>& b2a,    \
                                                                    const T              crustA, \
                                                                    const T              crustB, \
                                                                    Vector3<T>&          pa,     \
                                                                    Vector3<T>&          pb,     \
                                                                    uint&                nbIter);               \
    template __HOSTDEVICE__ T computeClosestPoints_GJK<T, GJK, ACC>(const Convex<T>&     a,      \
                                                                    const Convex<T>&     b,      \
                                                                    const Transform3<T>& a2w,    \
                                                                    const Transform3<T>& b2w,    \
                                                                    const T              crustA, \
                                                                    const T              crustB, \
                                                                    Vector3<T>&          pa,     \
                                                                    Vector3<T>&          pb,     \
                                                                    uint&                nbIter);               \
    template __HOSTDEVICE__ T computeClosestPoints_GJK<T, GJK, ACC>(const Convex<T>&     a,      \
                                                                    const Convex<T>&     b,      \
                                                                    const Vector3<T>&    v_b2a,  \
                                                                    const Quaternion<T>& q_b2a,  \
                                                                    const T              crustA, \
                                                                    const T              crustB, \
                                                                    Vector3<T>&          pa,     \
                                                                    Vector3<T>&          pb,     \
                                                                    uint&                nbIter);               \
    template __HOSTDEVICE__ T computeClosestPoints_GJK<T, GJK, ACC>(const Convex<T>&     a,      \
                                                                    const Convex<T>&     b,      \
                                                                    const Vector3<T>&    v_a2w,  \
                                                                    const Vector3<T>&    v_b2w,  \
                                                                    const Quaternion<T>& q_a2w,  \
                                                                    const Quaternion<T>& q_b2w,  \
                                                                    const T              crustA, \
                                                                    const T              crustB, \
                                                                    Vector3<T>&          pa,     \
                                                                    Vector3<T>&          pb,     \
                                                                    uint&                nbIter);               \
    template __HOSTDEVICE__ T computeClosestPoints_GJK<T, GJK, ACC>(const ShapeData<T>&  sdA,    \
                                                                    const ShapeData<T>&  sdB,    \
                                                                    const Vector3<T>&    v_b2a,  \
                                                                    const Quaternion<T>& q_b2a,  \
                                                                    const T              crustA, \
                                                                    const T              crustB, \
                                                                    Vector3<T>&          pa,     \
                                                                    Vector3<T>&          pb,     \
                                                                    uint&                nbIter);               \
    template __HOSTDEVICE__ T computeClosestPoints_GJK<T, GJK, ACC>(const ShapeData<T>&  sdA,    \
                                                                    const ShapeData<T>&  sdB,    \
                                                                    const Vector3<T>&    v_a2w,  \
                                                                    const Vector3<T>&    v_b2w,  \
                                                                    const Quaternion<T>& q_a2w,  \
                                                                    const Quaternion<T>& q_b2w,  \
                                                                    const T              crustA, \
                                                                    const T              crustB, \
                                                                    Vector3<T>&          pa,     \
                                                                    Vector3<T>&          pb,     \
                                                                    uint&                nbIter);
X(float, GJKType::JOHNSON, false)
X(double, GJKType::JOHNSON, false)
X(float, GJKType::JOHNSON, true)
X(double, GJKType::JOHNSON, true)
X(float, GJKType::SIGNEDVOLUME, false)
X(double, GJKType::SIGNEDVOLUME, false)
X(float, GJKType::SIGNEDVOLUME, true)
X(double, GJKType::SIGNEDVOLUME, true)
#undef X
