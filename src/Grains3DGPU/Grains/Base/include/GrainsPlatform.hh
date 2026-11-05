#ifndef GRAINS_PLATFORM_HH
#define GRAINS_PLATFORM_HH

#include "GrainsPlatform.hh"

// Compiler capability, not the selected runtime backend. Ordinary C++ translation units
// must not parse CUDA syntax or require CUDA headers, even in a CUDA-enabled build.
#if defined(__NVCC__) || defined(__CUDACC__)
#define GRAINS_CUDA_COMPILER 1
#else
#define GRAINS_CUDA_COMPILER 0
#endif

/** @brief Compiler macros */
#if GRAINS_CUDA_COMPILER
#define __HOST__ __host__
#define __DEVICE__ __device__
#define __HOSTDEVICE__ __host__ __device__
#define __MANAGED__ __managed__
#define __GLOBAL__ __global__
#define INLINE __inline__
#define __RESTRICT__ __restrict__
#else
// Host annotations for ordinary C++ compilation
#define __HOST__
#define __DEVICE__
#define __HOSTDEVICE__
#define __GLOBAL__
#define __MANAGED__
#define INLINE inline
#define __RESTRICT__
#endif

#endif
