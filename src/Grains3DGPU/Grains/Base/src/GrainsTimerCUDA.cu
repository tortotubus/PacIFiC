#include "GrainsTimerBackend.hh"

#include <cuda_runtime.h>

namespace grains::detail
{
    void synchronizeTimerDevice() noexcept
    {
        cudaDeviceSynchronize();
    }
}
