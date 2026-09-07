#ifndef _OBSTACLELOADING_HH_
#define _OBSTACLELOADING_HH_

#include "Vector3.hh"

// =================================================================================================
/** @brief Per-event record describing one prescribed-motion interval for one obstacle.

    For obstacle i, the host constructs one ObstacleMotionEvent entry per prescribed-motion interval

    @author A.Yazdani - 2026 - Construction */
// =================================================================================================
template <typename T>
struct ObstacleMotionEvent
{
    uint       obstacleId;       ///< Index into position/quaternion arrays
    T          tStart;           ///< Interval start time (inclusive)
    T          tEnd;             ///< Interval end time (exclusive)
    Vector3<T> linearVelocity;   ///< World-frame translational velocity [m/s]
    Vector3<T> angularVelocity;  ///< World-frame angular velocity [rad/s]
};

#endif
