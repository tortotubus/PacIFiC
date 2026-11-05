#ifndef _GRAINSUTILS_HH_
#define _GRAINSUTILS_HH_

#include "GrainsPlatform.hh"

#include "Basic.hh"
#include "GrainsAssert.hh"
#include "Vector3.hh"
#include <unistd.h>

#if GRAINS_CUDA_COMPILER
#include "GrainsCUDAUtils.hh"
#endif

// =================================================================================================
/** @brief Miscellaneous functionalities (mostly low-level) for Grains.

    @author A.Yazdani - 2025 - Construction */
// =================================================================================================
/** @name Miscellaneous functions and utilities for Grains */
//@{
// -------------------------------------------------------------------------------------------------
/** @brief Returns the available memory on the host in bytes */
__HOST__ static INLINE size_t getAvailableHostMemory()
{
    long pages     = sysconf(_SC_AVPHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
    return pages * page_size;
}

// -------------------------------------------------------------------------------------------------
/** @brief Writes a real number with a prescribed number of digits in a string
    @param figure the float number
    @param size number of digits */
template <typename T>
__HOST__ static INLINE std::string realToString(const T& figure, const int size)
{
    std::ostringstream oss;
    oss.width(size);
    oss << std::left << figure;
    return (oss.str());
}

// -------------------------------------------------------------------------------------------------
/** @brief Writes a float number with a prescribed format and a prescribed number of digits after
    the decimal point in a string.
    @param format the format
    @param digits number of digits after the decimal point
    @param number the float number */
template <typename T>
__HOST__ static INLINE std::string
    realToString(std::ios_base::fmtflags format, const int digits, const T& number)
{
    std::ostringstream oss;
    if(number != T(0))
    {
        oss.setf(format, std::ios::floatfield);
        oss.precision(digits);
    }
    oss << number;
    return (oss.str());
}

// -------------------------------------------------------------------------------------------------
/** @brief Writes a vector3 object in a string
    @param vec the vector3 object */
template <typename T>
__HOST__ static INLINE std::string Vector3ToString(const Vector3<T>& vec)
{
    std::ostringstream oss;
    oss << vec;
    return (oss.str());
}

// -------------------------------------------------------------------------------------------------
/** @brief Writes a message to stdout
    @param args the output messages */
template <typename... Args>
__HOST__ static constexpr INLINE void Gout(const Args&... args)
{
    ((std::cout << args << " "), ...);
    std::cout << std::endl;
}

template <typename... Args>
__HOST__ INLINE void GoutWI(const int numShift, const Args&... args)
{
    // indent
    auto shift = [](int n) { return std::string(n, ' '); };

    std::cout << shift(numShift);
    ((std::cout << args << " "), ...);
    std::cout << std::endl;
}

// -------------------------------------------------------------------------------------------------
/** @brief Computes a triangular hash for a pair of unsigned integers. Maps two IDs (x, y) to a
    unique single ID using upper-triangular indexing:
    hash = l * (l + 1) / 2 + s, where l = max(x,y) and s = min(x,y).
    This produces a symmetric hash (order-independent) for pairs.
    @param x first ID
    @param y second ID */
__HOSTDEVICE__ static INLINE uint triangularHash(uint x, uint y)
{
    uint s = min(x, y);  // smaller one
    uint l = max(x, y);  // larger one
    return (l * (l + 1) / 2 + s);
}

// -------------------------------------------------------------------------------------------------
/** @brief Computes a hash value for a pair of IDs using prime number multiplication. This provides
    good distribution of hash values to minimize collisions.
    @param k pair of IDs (i,j) where i < j */
__HOSTDEVICE__ static INLINE uint primeHash(IndexPair k)
{
    // Simple, fast hash using prime numbers for good distribution
    if(k.x > k.y)
    {
        const uint first = k.x;
        k.x = k.y;
        k.y = first;
    }
    uint h = k.x * 73856093u ^ k.y * 19349663u;
    return h;
}
//@}

#endif
