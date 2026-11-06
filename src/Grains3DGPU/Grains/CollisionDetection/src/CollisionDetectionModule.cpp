#include "GrainsUtils.hh"

#include "CollisionDetectionModuleImpl.hh"

namespace grains::cpu
{
    std::size_t getAvailableCollisionMemory()
    {
        return getAvailableHostMemory();
    }
}

template class CollisionDetectionModule<float, MemType::HOST>;
template class CollisionDetectionModule<double, MemType::HOST>;
