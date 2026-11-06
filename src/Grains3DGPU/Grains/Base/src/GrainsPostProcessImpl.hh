#ifndef GRAINS_POST_PROCESS_IMPL_HH
#define GRAINS_POST_PROCESS_IMPL_HH

#include "Grains.hh"

template <typename T>
template <MemType M>
void Grains<T>::postProcess(const std::unique_ptr<ComponentManager<T, M>>& cm)
{
    using GP = GrainsParameters<T>;
    auto& state = GP::m_simulationState;

    if(GP::m_tSave.empty())
        return;

    if(GP::m_tSave.front() - state.time < T(0.01) * GP::m_dt)
    {
        GP::m_tSave.pop();
        if constexpr(M == MemType::DEVICE)
        {
            cm->copyTo_PostProcessing(m_components);
            for(auto& writer : m_postProcessor)
                writer->PostProcessing(m_rigidBodyList, m_components, state.time);
        }
        else
        {
            for(auto& writer : m_postProcessor)
                writer->PostProcessing(m_rigidBodyList, cm, state.time);
        }
    }

    if(!GP::m_tSave.empty() && state.time > GP::m_tSave.front())
        GP::m_tSave.pop();
}

#endif
