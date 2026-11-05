#include "ForceModule.hh"
#include "ForceComputation.hh"
#include "ForceComputationCPU.hh"
#include "ForceComputationCUDA.hh"

// CUDA resources are created/destroyed by the CUDA backend, never by host modules.
template <typename T, MemType M>
ForceModule<T, M>::ForceModule(size_t pairCapacity, bool isContactWithMemory)
{
    if constexpr(M == MemType::DEVICE)
        m_backend.state = grains::cuda::createForceState<T>(pairCapacity);
    if(isContactWithMemory)
    {
        uint hashCapacity = static_cast<uint>(pairCapacity / 0.7);
        m_contactTable.allocate(hashCapacity, static_cast<uint>(pairCapacity));
    }
}

template <typename T, MemType M>
void ForceModule<T, M>::resizeBuffers(size_t newPairCapacity)
{
    if constexpr(M == MemType::DEVICE)
        grains::cuda::resizeForceState(*m_backend.state, newPairCapacity);
}

// -------------------------------------------------------------------------------------------------
// Periodic mark-and-sweep
template <typename T, MemType M>
void ForceModule<T, M>::cleanupContactTable()
{
    using GP = GrainsParameters<T>;
    if(!GP::m_isContactWithMemory)
        return;

    auto& SS = GP::m_simulationState;
    if(SS.neighborListUpdateCount % 1000 == 0)
        m_contactTable.markAndSweep();
}

// Accumulate dependent sub-body torces into their composite masters.
template <typename T, MemType M>
void ForceModule<T, M>::assembleCompositeTorces(GrainsMemBuffer<Torce<T>, M>&         torce,
                                                const GrainsMemBuffer<Vector3<T>, M>& position,
                                                const GrainsMemBuffer<uint, M>&       bodyTag,
                                                const GrainsMemBuffer<uint, M>&       masterSlot,
                                                const ComponentCounts&                counts)
{
    if(counts.numSubBodies == 0)
        return;
    const CompositeForceView<T> view{torce.getData(),
                                     position.getData(),
                                     masterSlot.getData(),
                                     bodyTag.getData(),
                                     M == MemType::HOST ? counts.numObstacles : 0u,
                                     counts.numObstacles + counts.numParticles};
    if constexpr(M == MemType::HOST)
        grains::cpu::assembleCompositeForces(view);
    else
        grains::cuda::assembleCompositeForces(view);
}

template <typename T, MemType M>
void ForceModule<T, M>::run(const GrainsMemBuffer<ContactForceModel<T>*, M>& CF,
                            const GrainsMemBuffer<RigidBody<T>*, M>*         rigidBody,
                            const GrainsMemBuffer<Vector3<T>, M>&            position,
                            const GrainsMemBuffer<Kinematics<T>, M>&         velocity,
                            const GrainsMemBuffer<IndexPair, M>&             pairList,
                            const GrainsMemBuffer<ContactInfo<T>, M>&        contactInfo,
                            GrainsMemBuffer<Torce<T>, M>&                    torce,
                            const GrainsMemBuffer<uint, M>&                  bodyTag,
                            const GrainsMemBuffer<uint, M>&                  masterSlot,
                            const ComponentCounts&                           counts)
{
    using GP    = GrainsParameters<T>;
    auto& timer = GP::m_fmTimer;
    timer.start(FMStage::Total);
    cleanupContactTable();

    const ContactForceView<T> contacts{CF.getData(),
                                       pairList.getData(),
                                       contactInfo.getData(),
                                       position.getData(),
                                       velocity.getData(),
                                       torce.getData(),
                                       m_contactTable.getView(),
                                       counts.numPairs};
    if constexpr(M == MemType::HOST)
    {
        timer.start(FMStage::ComputeForces);
        grains::cpu::computeContactForces(contacts);
        timer.stop(FMStage::ComputeForces);
    }
    else
        grains::cuda::computeContactForces(*m_backend.state, contacts, GP::m_useCompaction, timer);

    const ExternalForceView<T> external{rigidBody->getData(),
                                        torce.getData(),
                                        GP::m_gravity[X],
                                        GP::m_gravity[Y],
                                        GP::m_gravity[Z],
                                        counts.numObstacles,
                                        counts.numObstacles + counts.numParticles};
    timer.start(FMStage::ExternalForces);
    if constexpr(M == MemType::HOST)
        grains::cpu::addExternalForces(external);
    else
        grains::cuda::addExternalForces(external);
    timer.stop(FMStage::ExternalForces);

    if(counts.numSubBodies > 0)
    {
        timer.start(FMStage::AssembleComposites);
        assembleCompositeTorces(torce, position, bodyTag, masterSlot, counts);
        timer.stop(FMStage::AssembleComposites);
    }
    timer.stop(FMStage::Total);
}

template class ForceModule<float, MemType::HOST>;
template class ForceModule<double, MemType::HOST>;
#if GRAINS_CUDA_COMPILER
template class ForceModule<float, MemType::DEVICE>;
template class ForceModule<double, MemType::DEVICE>;
#endif
