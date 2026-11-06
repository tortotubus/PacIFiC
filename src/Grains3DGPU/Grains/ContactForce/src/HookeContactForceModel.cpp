#include "HookeContactForceModelImpl.hh"

// -------------------------------------------------------------------------------------------------
// Constructor with an XML node
template <typename T>
__HOST__ HookeContactForceModel<T>::HookeContactForceModel(grains::xml::OptionalElement root)
{
    const auto& element = grains::xml::requireElement(root);
    m_kn   = T(element.requireAttribute<double>("kn"));
    m_en   = T(element.requireAttribute<double>("en"));
    m_etat = T(element.requireAttribute<double>("etat"));
    m_muc  = T(element.requireAttribute<double>("muc"));
    m_kr   = T(element.requireAttribute<double>("kr"));

    m_muen = log(m_en) / sqrt(PI<T> * PI<T> + log(m_en) * log(m_en));
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class HookeContactForceModel<float>;
template class HookeContactForceModel<double>;
