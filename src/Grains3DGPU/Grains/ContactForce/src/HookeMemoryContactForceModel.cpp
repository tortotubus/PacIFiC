#include "HookeMemoryContactForceModelImpl.hh"

#include "GrainsParameters.hh"

// -------------------------------------------------------------------------------------------------
// Constructor with an XML node
template <typename T>
__HOST__ HookeMemoryContactForceModel<T>::HookeMemoryContactForceModel(
    grains::xml::OptionalElement root)
{
    GAssert(GrainsParameters<T>::m_dt > T(0), "Time step not defined! Aborting Grains!");

    const auto& element = grains::xml::requireElement(root);
    m_kn   = T(element.requireAttribute<double>("kn"));
    m_en   = T(element.requireAttribute<double>("en"));
    m_kt   = T(element.requireAttribute<double>("kt"));
    m_etat = T(element.requireAttribute<double>("etat"));
    m_muc  = T(element.requireAttribute<double>("muc"));
    if(auto mur = element.attributeAs<double>("mur"))
        m_mur = T(*mur);
    if(auto etarpf = element.attributeAs<double>("etarpf"))
        m_etarpf = T(*etarpf);

    m_muen = log(m_en) / sqrt(PI<T> * PI<T> + log(m_en) * log(m_en));
    m_dt   = GrainsParameters<T>::m_dt;
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class HookeMemoryContactForceModel<float>;
template class HookeMemoryContactForceModel<double>;
