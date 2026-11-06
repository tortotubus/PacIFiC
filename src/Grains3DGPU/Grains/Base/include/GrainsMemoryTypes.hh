#ifndef GRAINS_MEMORY_TYPES_HH
#define GRAINS_MEMORY_TYPES_HH

#include <cstddef>
#include <cstdint>

enum class MemType
{
    HOST,
    DEVICE,
    MANAGED,
    PINNED,
    MAPPED,
    UNKNOWN
};

/** Backend-neutral native stream handle. CUDA backends reinterpret it as cudaStream_t. */
using GrainsStream = void*;

enum class GrainsScalarKind : std::uint8_t
{
    SIGNED_INTEGER,
    UNSIGNED_INTEGER,
    FLOATING_POINT,
    BOOLEAN
};

/** Small trivially-copyable scalar packet passed across the C++/CUDA backend boundary. */
struct alignas(16) GrainsScalarValue
{
    unsigned char    bytes[16]{};
    std::uint8_t     size = 0;
    GrainsScalarKind kind = GrainsScalarKind::UNSIGNED_INTEGER;
};

#endif
