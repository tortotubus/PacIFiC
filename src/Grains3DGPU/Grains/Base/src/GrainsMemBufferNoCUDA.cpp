#include "GrainsMemBufferBackend.hh"

#include <cstdlib>

#include "GrainsAssert.hh"

namespace
{
    [[noreturn]] void cudaRequired()
    {
        GAbort("GrainsMemBuffer operation requires a CUDA-enabled build");
        std::abort();
    }
}

namespace grains::memory
{
    Allocation allocate(MemType, std::size_t)
    {
        cudaRequired();
    }

    void release(MemType, Allocation allocation)
    {
        if(allocation.pointer == nullptr)
            return;
        cudaRequired();
    }

    void copyAsync(void*, MemType, const void*, MemType, std::size_t, GrainsStream)
    {
        cudaRequired();
    }

    void synchronize(GrainsStream)
    {
        cudaRequired();
    }

    void setBytesAsync(void*, int, std::size_t, GrainsStream)
    {
        cudaRequired();
    }

    void fillScalarAsync(void*, std::size_t, const GrainsScalarValue&, GrainsStream)
    {
        cudaRequired();
    }

    void sequenceScalarAsync(void*, std::size_t, const GrainsScalarValue&, GrainsStream)
    {
        cudaRequired();
    }

    void freeDevicePointers(void*, std::size_t)
    {
        cudaRequired();
    }
}
