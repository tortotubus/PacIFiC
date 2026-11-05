#ifndef GRAINS_TYPES_HH
#define GRAINS_TYPES_HH

#include "GrainsPlatform.hh"
#if GRAINS_CUDA_COMPILER
#include <vector_types.h>
#endif

// Simulation indices keep the CUDA layout when compiled by CUDA, and use plain
// data structures on the host. No CUDA runtime types enter a C++ translation unit.
using uint = unsigned int;

#if GRAINS_CUDA_COMPILER
using IndexPair      = uint2;
using CellIndex      = uint3;
using CellDimensions = uint4;
#else
struct CellIndex
{
    uint x, y, z;
};
struct alignas(16) CellDimensions
{
    uint x, y, z, w;
};
struct alignas(8) IndexPair
{
    uint x;
    uint y;
};
#endif

static_assert(sizeof(IndexPair) == 8 && alignof(IndexPair) == 8);
static_assert(sizeof(CellIndex) == 12 && alignof(CellIndex) == 4);
static_assert(sizeof(CellDimensions) == 16 && alignof(CellDimensions) == 16);

__HOSTDEVICE__ inline CellIndex makeCellIndex(uint x, uint y, uint z)
{
    return {x, y, z};
}

__HOSTDEVICE__ inline CellDimensions makeCellDimensions(uint x, uint y, uint z, uint w)
{
    return {x, y, z, w};
}

__HOSTDEVICE__ inline IndexPair makeIndexPair(uint x, uint y)
{
    return {x, y};
}
#endif
