#include "GrainsMemBufferBackend.hh"

#include <cuda_runtime.h>

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <new>
#include <vector>

#include "GrainsAssert.hh"
#include "GrainsCUDAError.hh"

namespace
{
    cudaStream_t nativeStream(GrainsStream stream)
    {
        return reinterpret_cast<cudaStream_t>(stream);
    }

    cudaMemcpyKind copyKind(MemType source, MemType destination)
    {
        if(source == MemType::HOST && destination == MemType::HOST)
            return cudaMemcpyHostToHost;
        if(source == MemType::HOST && destination == MemType::DEVICE)
            return cudaMemcpyHostToDevice;
        if(source == MemType::DEVICE && destination == MemType::HOST)
            return cudaMemcpyDeviceToHost;
        if(source == MemType::DEVICE && destination == MemType::DEVICE)
            return cudaMemcpyDeviceToDevice;
        if(source == MemType::HOST && destination == MemType::MANAGED)
            return cudaMemcpyHostToDevice;
        if(source == MemType::MANAGED && destination == MemType::HOST)
            return cudaMemcpyDeviceToHost;
        if(source == MemType::DEVICE && destination == MemType::MANAGED)
            return cudaMemcpyDeviceToDevice;
        if(source == MemType::MANAGED && destination == MemType::DEVICE)
            return cudaMemcpyDeviceToDevice;
        if(source == MemType::MANAGED && destination == MemType::MANAGED)
            return cudaMemcpyDeviceToDevice;

        const bool sourceHostMemory
            = source == MemType::HOST || source == MemType::PINNED || source == MemType::MAPPED;
        const bool destinationHostMemory = destination == MemType::HOST
                                           || destination == MemType::PINNED
                                           || destination == MemType::MAPPED;
        if(sourceHostMemory && destinationHostMemory)
            return cudaMemcpyHostToHost;
        if(sourceHostMemory && destination == MemType::DEVICE)
            return cudaMemcpyHostToDevice;
        if(source == MemType::DEVICE && destinationHostMemory)
            return cudaMemcpyDeviceToHost;
        return cudaMemcpyDefault;
    }

    template <typename Bits>
    __GLOBAL__ void fillBitsKernel(Bits* destination, std::size_t count, Bits value)
    {
        const std::size_t index = blockIdx.x * blockDim.x + threadIdx.x;
        if(index < count)
            destination[index] = value;
    }

    struct Bits128
    {
        std::uint64_t low;
        std::uint64_t high;
    };

    template <typename T>
    T unpack(const GrainsScalarValue& value)
    {
        T result{};
        std::memcpy(&result, value.bytes, sizeof(T));
        return result;
    }

    template <typename T>
    void launchFill(void*                    destination,
                    std::size_t              count,
                    const GrainsScalarValue& value,
                    cudaStream_t             stream)
    {
        constexpr unsigned int threads = 256;
        const unsigned int     blocks  = static_cast<unsigned int>((count + threads - 1) / threads);
        fillBitsKernel<<<blocks, threads, 0, stream>>>(static_cast<T*>(destination),
                                                       count,
                                                       unpack<T>(value));
        cudaErrCheck(cudaGetLastError());
    }

    template <typename T>
    __GLOBAL__ void sequenceKernel(T* destination, std::size_t count, T start)
    {
        const std::size_t index = blockIdx.x * blockDim.x + threadIdx.x;
        if(index < count)
            destination[index] = static_cast<T>(start + index);
    }

    template <typename T>
    void launchSequence(void*                    destination,
                        std::size_t              count,
                        const GrainsScalarValue& start,
                        cudaStream_t             stream)
    {
        constexpr unsigned int threads = 256;
        const unsigned int     blocks  = static_cast<unsigned int>((count + threads - 1) / threads);
        sequenceKernel<<<blocks, threads, 0, stream>>>(static_cast<T*>(destination),
                                                       count,
                                                       unpack<T>(start));
        cudaErrCheck(cudaGetLastError());
    }

    __GLOBAL__ void sequenceBoolKernel(bool* destination, std::size_t count, bool start)
    {
        const std::size_t index = blockIdx.x * blockDim.x + threadIdx.x;
        if(index < count)
            destination[index] = static_cast<bool>(start + index);
    }
}

namespace grains::memory
{
    Allocation allocate(MemType type, std::size_t bytes)
    {
        Allocation allocation;
        switch(type)
        {
        case MemType::HOST:
            allocation.pointer = std::malloc(bytes);
            if(allocation.pointer == nullptr)
                throw std::bad_alloc();
            break;
        case MemType::DEVICE:
            cudaErrCheck(cudaMalloc(&allocation.pointer, bytes));
            break;
        case MemType::MANAGED:
            cudaErrCheck(cudaMallocManaged(&allocation.pointer, bytes));
            break;
        case MemType::PINNED:
            cudaErrCheck(cudaMallocHost(&allocation.pointer, bytes));
            break;
        case MemType::MAPPED:
            cudaErrCheck(cudaHostAlloc(&allocation.pointer, bytes, cudaHostAllocMapped));
            cudaErrCheck(cudaHostGetDevicePointer(&allocation.deviceAlias, allocation.pointer, 0));
            break;
        case MemType::UNKNOWN:
            GAbort("Cannot allocate UNKNOWN memory type");
        }
        return allocation;
    }

    void release(MemType type, Allocation allocation)
    {
        if(allocation.pointer == nullptr)
            return;

        switch(type)
        {
        case MemType::HOST:
            std::free(allocation.pointer);
            break;
        case MemType::DEVICE:
        case MemType::MANAGED:
            cudaErrCheck(cudaFree(allocation.pointer));
            break;
        case MemType::PINNED:
        case MemType::MAPPED:
            cudaErrCheck(cudaFreeHost(allocation.pointer));
            break;
        case MemType::UNKNOWN:
            GAbort("Cannot release UNKNOWN memory type");
        }
    }

    void copyAsync(void*        destination,
                   MemType      destinationType,
                   const void*  source,
                   MemType      sourceType,
                   std::size_t  bytes,
                   GrainsStream stream)
    {
        cudaErrCheck(cudaMemcpyAsync(destination,
                                     source,
                                     bytes,
                                     copyKind(sourceType, destinationType),
                                     nativeStream(stream)));
    }

    void synchronize(GrainsStream stream)
    {
        cudaErrCheck(cudaStreamSynchronize(nativeStream(stream)));
    }

    void setBytesAsync(void* destination, int value, std::size_t bytes, GrainsStream stream)
    {
        cudaErrCheck(cudaMemsetAsync(destination, value, bytes, nativeStream(stream)));
    }

    void fillScalarAsync(void*                    destination,
                         std::size_t              count,
                         const GrainsScalarValue& value,
                         GrainsStream             stream)
    {
        switch(value.size)
        {
        case 1:
            launchFill<std::uint8_t>(destination, count, value, nativeStream(stream));
            break;
        case 2:
            launchFill<std::uint16_t>(destination, count, value, nativeStream(stream));
            break;
        case 4:
            launchFill<std::uint32_t>(destination, count, value, nativeStream(stream));
            break;
        case 8:
            launchFill<std::uint64_t>(destination, count, value, nativeStream(stream));
            break;
        case 16:
            launchFill<Bits128>(destination, count, value, nativeStream(stream));
            break;
        default:
            GAbort("Unsupported scalar byte width for GrainsMemBuffer::fill:", value.size);
        }
    }

    void sequenceScalarAsync(void*                    destination,
                             std::size_t              count,
                             const GrainsScalarValue& start,
                             GrainsStream             stream)
    {
        const cudaStream_t cudaStream = nativeStream(stream);
        if(start.kind == GrainsScalarKind::BOOLEAN)
        {
            constexpr unsigned int threads = 256;
            const unsigned int blocks = static_cast<unsigned int>((count + threads - 1) / threads);
            sequenceBoolKernel<<<blocks, threads, 0, cudaStream>>>(static_cast<bool*>(destination),
                                                                   count,
                                                                   unpack<bool>(start));
            cudaErrCheck(cudaGetLastError());
            return;
        }

        if(start.kind == GrainsScalarKind::FLOATING_POINT)
        {
            if(start.size == sizeof(float))
                launchSequence<float>(destination, count, start, cudaStream);
            else if(start.size == sizeof(double))
                launchSequence<double>(destination, count, start, cudaStream);
            else
                GAbort("Unsupported floating-point width for GrainsMemBuffer::sequence:",
                       start.size);
            return;
        }

        const bool isSigned = start.kind == GrainsScalarKind::SIGNED_INTEGER;
        switch(start.size)
        {
        case 1:
            isSigned ? launchSequence<std::int8_t>(destination, count, start, cudaStream)
                     : launchSequence<std::uint8_t>(destination, count, start, cudaStream);
            break;
        case 2:
            isSigned ? launchSequence<std::int16_t>(destination, count, start, cudaStream)
                     : launchSequence<std::uint16_t>(destination, count, start, cudaStream);
            break;
        case 4:
            isSigned ? launchSequence<std::int32_t>(destination, count, start, cudaStream)
                     : launchSequence<std::uint32_t>(destination, count, start, cudaStream);
            break;
        case 8:
            isSigned ? launchSequence<std::int64_t>(destination, count, start, cudaStream)
                     : launchSequence<std::uint64_t>(destination, count, start, cudaStream);
            break;
        default:
            GAbort("Unsupported integer width for GrainsMemBuffer::sequence:", start.size);
        }
    }

    void freeDevicePointers(void* pointerArray, std::size_t count)
    {
        std::vector<void*> hostPointers(count);
        cudaErrCheck(cudaMemcpy(hostPointers.data(),
                                pointerArray,
                                count * sizeof(void*),
                                cudaMemcpyDeviceToHost));
        for(void* pointer : hostPointers)
        {
            if(pointer != nullptr)
                cudaErrCheck(cudaFree(pointer));
        }
        cudaErrCheck(cudaMemset(pointerArray, 0, count * sizeof(void*)));
        cudaErrCheck(cudaDeviceSynchronize());
    }
}
