#ifndef GRAINS_ASSERT_HH
#define GRAINS_ASSERT_HH

#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

#include "GrainsPlatform.hh"

// Helper functions for device printf with different types
__DEVICE__ INLINE void print_device_arg(const char* arg)
{
    printf("%s ", arg);
}
__DEVICE__ INLINE void print_device_arg(char* arg)
{
    printf("%s ", arg);
}
__DEVICE__ INLINE void print_device_arg(size_t arg)
{
    printf("%zu ", arg);
}
__DEVICE__ INLINE void print_device_arg(int arg)
{
    printf("%d ", arg);
}
__DEVICE__ INLINE void print_device_arg(unsigned int arg)
{
    printf("%u ", arg);
}
__DEVICE__ INLINE void print_device_arg(long arg)
{
    printf("%ld ", arg);
}
__DEVICE__ INLINE void print_device_arg(float arg)
{
    printf("%f ", arg);
}
__DEVICE__ INLINE void print_device_arg(double arg)
{
    printf("%f ", arg);
}
__HOST__ INLINE void print_device_arg(const std::string& arg)
{
    printf("%s ", arg.c_str());
}

template <typename... Args>
__HOSTDEVICE__ INLINE void GAbort(const Args&... args)
{
#ifdef __CUDA_ARCH__
    printf("[DEVICE] ");
    (print_device_arg(args), ...);
    printf("\n");
    __trap();  // aborts the kernel
#else
    std::cerr << "[HOST] ";
    ((std::cerr << args << " "), ...);
    std::cerr << std::endl;
    std::abort();
#endif
}

// -------------------------------------------------------------------------------------------------
/** @brief Assert function that aborts the program if the condition is false
    @param condition the condition to check
    @param args the message(s) to display if the assertion fails */
template <typename... Args>
__HOSTDEVICE__ INLINE void GAssert(bool condition, const Args&... args)
{
    if(!condition)
        GAbort("GAssert failed:", args...);
}

#endif
