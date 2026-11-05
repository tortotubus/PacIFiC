#ifndef GRAINS_CUDA_ERROR_HH
#define GRAINS_CUDA_ERROR_HH

#include "GrainsPlatform.hh"
#if !GRAINS_CUDA_COMPILER
#error "GrainsCUDAError.hh requires a CUDA compiler"
#endif
#include <cstdio>
#include <cstdlib>
#include <cuda_runtime.h>

// Keep runtime error reporting usable without the legacy Basic.hh namespace imports.
// Macro for outputting CUDA errors
#define cudaErrCheck(ans) cudaAssert((ans), __FILE__, __LINE__);

/** @brief Returns CUDA error.
    @param code the error code
    @param file the file name
    @param line the line number
    @param abort whether to abort the program */
__HOST__ static INLINE void
    cudaAssert(cudaError_t code, const char* file, int line, bool abort = false)
{
    if(code != cudaSuccess)
    {
        fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
        if(abort)
            exit(code);
    }
}

#endif
