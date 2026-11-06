#ifndef GRAINS_TYPES_HH
#define GRAINS_TYPES_HH

#include "GrainsPlatform.hh"

// Project-owned index types keep one name and ABI across C++ and CUDA translation units.
using uint = unsigned int;

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
