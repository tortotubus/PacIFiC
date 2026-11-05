#ifndef _PARTICLESORTER_HH_
#define _PARTICLESORTER_HH_

#include "GrainsMemBuffer.hh"
#include "GrainsParameters.hh"
#include "ParticleSorting.hh"
#include "ParticleSortingCPU.hh"
#include "ParticleSortingCUDA.hh"

namespace grains::detail
{
    template <typename T, MemType M>
    struct ParticleSortStorage;
    template <typename T>
    struct ParticleSortStorage<T, MemType::HOST>
    {
        cpu::ParticleSortStatePtr<T> state;
    };
    template <typename T>
    struct ParticleSortStorage<T, MemType::DEVICE>
    {
        cuda::ParticleSortStatePtr<T> state;
    };
}

/** @brief Reorders the particle range by Morton cell codes, leaving obstacles in place.
    Backends own their cells and reusable sort/gather storage. CUDA resources and
    launches are implemented in .cu sources; this wrapper is ordinary C++. */
template <typename T, MemType M>
class ParticleSorter
{
private:
    grains::detail::ParticleSortStorage<T, M> m_backend;

public:
    ParticleSorter() = default;
    ParticleSorter(uint numObstacles, uint numParticles)
    {
        const auto& params = GrainsParameters<T>::m_collisionDetection.linkedCellParameters;
        if constexpr(M == MemType::HOST)
            m_backend.state = grains::cpu::createParticleSortState<T>(params.minCorner,
                                                                      params.maxCorner,
                                                                      params.cellSizeFactor,
                                                                      numParticles);
        else
            m_backend.state = grains::cuda::createParticleSortState<T>(params.minCorner,
                                                                       params.maxCorner,
                                                                       params.cellSizeFactor,
                                                                       numParticles);
    }
    ~ParticleSorter()                                = default;
    ParticleSorter(const ParticleSorter&)            = delete;
    ParticleSorter& operator=(const ParticleSorter&) = delete;
    ParticleSorter(ParticleSorter&&)                 = default;
    ParticleSorter& operator=(ParticleSorter&&)      = default;

    void sortParticles(GrainsMemBuffer<Vector3<T>, M>&    position,
                       GrainsMemBuffer<Kinematics<T>, M>& velocity,
                       GrainsMemBuffer<Quaternion<T>, M>& quaternion,
                       GrainsMemBuffer<Torce<T>, M>&      torce,
                       GrainsMemBuffer<uint, M>&          bodyTag,
                       GrainsMemBuffer<Vector3<T>, M>&    localPos,
                       GrainsMemBuffer<Quaternion<T>, M>& localQuat,
                       uint                               numObstacles,
                       uint                               numParticles)
    {
        if(numParticles == 0)
            return;
        GAssert(m_backend.state != nullptr, "ParticleSorter must be initialized before sorting");
        const ParticleSortView<T> view{position.getData() + numObstacles,
                                       velocity.getData() + numObstacles,
                                       quaternion.getData() + numObstacles,
                                       torce.getData() + numObstacles,
                                       bodyTag.getData() + numObstacles,
                                       localPos.getData() + numObstacles,
                                       localQuat.getData() + numObstacles,
                                       numParticles};
        if constexpr(M == MemType::HOST)
            grains::cpu::sortParticles(*m_backend.state, view);
        else
            grains::cuda::sortParticles(*m_backend.state, view);
    }
};
#endif
