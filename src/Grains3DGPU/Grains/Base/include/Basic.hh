#ifndef _BASIC_HH_
#define _BASIC_HH_

#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include <algorithm>
#include <bits/stdc++.h>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <execution>
#include <float.h>
#include <ios>
#include <iostream>
#include <istream>
#include <limits>
#include <list>
#include <new>
#include <ostream>
#include <random>
#include <stdio.h>
#include <string>

#include "GrainsPlatform.hh"

#if GRAINS_CUDA_COMPILER
#include <cuda.h>
#include <cuda_runtime.h>
#include <curand.h>
#include <curand_kernel.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Xml.hh"

// The existing Grains3DGPU headers use unqualified standard-library names.
// ReaderXML.hh historically supplied this directive; retain that source
// compatibility here until those headers are cleaned up independently.
using namespace std;

// =================================================================================================
/** @brief Various constants and type definitions.

    @author A.Yazdani - 2023 - Construction */
// =================================================================================================
/** @name Macros */
//@{
// -------------------------------------------------------------------------------------------------
/** @brief Type macros */
//@}

/** @name Typedefs */
//@{
// -------------------------------------------------------------------------------------------------
/** @brief Unsigned integer type */
#include "GrainsTypes.hh"
//@}

/** @name Enumerations */
//@{
// -------------------------------------------------------------------------------------------------
/** @brief Space dimensions */
enum Direction
{
    X,    // x direction
    Y,    // y direction
    Z,    // z direction
    W,    // scalar component of quaternions
    NONE  // no direction
};

// -------------------------------------------------------------------------------------------------
/** @brief Matrix dimensions */
enum MatDirection
{
    XX,
    XY,
    XZ,
    YX,
    YY,
    YZ,
    ZX,
    ZY,
    ZZ
};
//@}

/** @name Constant expressions at compile-time */
//@{
// -------------------------------------------------------------------------------------------------
/** \brief PI value */
template <class T>
inline constexpr T PI = T(3.1415926535897932385L);

// -------------------------------------------------------------------------------------------------
/** \brief 2*PI value */
template <class T>
inline constexpr T TWO_PI = T(6.28318530717958623200L);

// -------------------------------------------------------------------------------------------------
/** \brief Degree per radian value */
template <class T>
inline constexpr T DEGS_PER_RAD = T(57.29577951308232286465L);

// -------------------------------------------------------------------------------------------------
/** \brief Radian per degree value */
template <class T>
inline constexpr T RADS_PER_DEG = T(0.01745329251994329547L);

// -------------------------------------------------------------------------------------------------
/** \brief High (Machine) epsilon value */
template <class T>
inline constexpr T HIGHEPS = T(1.e-15);
template <>
inline constexpr float HIGHEPS<float> = float(1.e-08);

// -------------------------------------------------------------------------------------------------
/** \brief Epsilon value */
template <class T>
inline constexpr T EPS = T(1.e-10);
template <>
inline constexpr float EPS<float> = float(1.e-05);

// -------------------------------------------------------------------------------------------------
/** \brief Low epsilon value */
template <class T>
inline constexpr T LOWEPS = T(1.e-05);
template <>
inline constexpr float LOWEPS<float> = float(1.e-03);
//@}

#endif
