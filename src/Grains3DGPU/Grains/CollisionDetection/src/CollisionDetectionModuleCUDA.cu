#include "GrainsCUDAUtils.hh"

#include "CollisionDetectionModuleImpl.hh"

namespace grains::cuda
{
    std::size_t getAvailableCollisionMemory()
    {
        return getAvailableDeviceMemory();
    }
}

template class CollisionDetectionModule<float, MemType::DEVICE>;
template class CollisionDetectionModule<double, MemType::DEVICE>;
