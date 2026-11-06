#ifndef GRAINS_CUDA_UTILS_HH
#define GRAINS_CUDA_UTILS_HH

#include "Basic.hh"
#include "GrainsParameters.hh"

#if !GRAINS_CUDA_COMPILER
#error "GrainsCUDAUtils.hh requires a CUDA compiler"
#endif

#include "GrainsCUDAError.hh"

// -------------------------------------------------------------------------------------------------
/** @brief Convert CUDA runtime properties to the backend-neutral launch-sizing representation. */
__HOST__ static INLINE GrainsDeviceProperties
    makeGrainsDeviceProperties(const cudaDeviceProp& properties)
{
    return {properties.multiProcessorCount, properties.maxThreadsPerBlock};
}

// -------------------------------------------------------------------------------------------------
/** @brief Query the launch-sizing properties for a CUDA device. */
__HOST__ static INLINE GrainsDeviceProperties getGrainsDeviceProperties(const int device = 0)
{
    cudaDeviceProp properties{};
    cudaErrCheck(cudaGetDeviceProperties(&properties, device));
    return makeGrainsDeviceProperties(properties);
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns the available memory on the device in bytes */
__HOST__ static INLINE size_t getAvailableDeviceMemory()
{
    size_t free_byte;
    size_t total_byte;
    cudaErrCheck(cudaMemGetInfo(&free_byte, &total_byte));
    return free_byte;
}

// -------------------------------------------------------------------------------------------------
/** @brief computes the optimal number of threads and blocks for a given number of elements and an
    architecture
    @param numElements the number of elements
    @param prop the device properties
    @param numBlocks the number of blocks (output)
    @param numThreads the number of threads per block (output) */
__HOST__ static INLINE void computeOptimalThreadsAndBlocks(const uint            numElements,
                                                           const GrainsDeviceProperties& prop,
                                                           uint&                 numBlocks,
                                                           uint&                 numThreads)
{
    if(numElements == 0)
    {
        fprintf(stderr, "[FATAL] computeOptimalThreadsAndBlocks called with numElements=0\n");
        abort();
    }
    constexpr uint maxThreads = 256;  // Avoid 1024 unless necessary
    constexpr uint minThreads = 32;
    const uint     minBlocks  = 2 * prop.multiProcessorCount;

    // Start with 128 threads and compute how many blocks we need
    numThreads = 128;
    numBlocks  = (numElements + numThreads - 1) / numThreads;

    // If we're not filling the SMs enough, reduce threads but ensure full coverage
    if(numBlocks < minBlocks && numThreads > minThreads)
    {
        numThreads = minThreads;
        numBlocks  = max(minBlocks, (numElements + numThreads - 1) / numThreads);
    }
    // If we're oversubscribing the SMs too much, increase thread count
    if(numBlocks > 4 * prop.multiProcessorCount && numThreads < maxThreads)
    {
        numThreads = maxThreads;
        numBlocks  = (numElements + numThreads - 1) / numThreads;
    }
}

// -------------------------------------------------------------------------------------------------
/** @brief Writes a IndexPair object in a string
    @param os the output stream
    @param v the IndexPair object */
__HOST__ static INLINE std::ostream& operator<<(std::ostream& os, const IndexPair& v)
{
    os << "(" << v.x << ", " << v.y << ")";
    return os;
}

#endif
