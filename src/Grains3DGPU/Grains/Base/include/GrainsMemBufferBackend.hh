#ifndef GRAINS_MEM_BUFFER_BACKEND_HH
#define GRAINS_MEM_BUFFER_BACKEND_HH

#include <cstddef>

#include "GrainsMemoryTypes.hh"

namespace grains::memory
{
    struct Allocation
    {
        void* pointer     = nullptr;
        void* deviceAlias = nullptr;
    };

    Allocation allocate(MemType type, std::size_t bytes);
    void       release(MemType type, Allocation allocation);

    void copyAsync(void*        destination,
                   MemType      destinationType,
                   const void*  source,
                   MemType      sourceType,
                   std::size_t  bytes,
                   GrainsStream stream);
    void synchronize(GrainsStream stream);

    void setBytesAsync(void* destination, int value, std::size_t bytes, GrainsStream stream);
    void fillScalarAsync(void*                    destination,
                         std::size_t              count,
                         const GrainsScalarValue& value,
                         GrainsStream             stream);
    void sequenceScalarAsync(void*                    destination,
                             std::size_t              count,
                             const GrainsScalarValue& start,
                             GrainsStream             stream);

    /** Preserve the legacy behavior: copy device pointers to the host and cudaFree each one. */
    void freeDevicePointers(void* pointerArray, std::size_t count);
}

#endif
