#include "RigidBodyImpl.hh"

// -------------------------------------------------------------------------------------------------
// Constructor with an XML input
template <typename T>
__HOST__ RigidBody<T>::RigidBody(grains::xml::OptionalElement root)
{
    const auto& element = grains::xml::requireElement(root);
    // Convex
    auto shape = element.getChild("Convex");
    m_convex       = ConvexFactory<T>::create(shape);
    // Crust thickness
    m_crustThickness = T(grains::xml::requireElement(shape)
                             .requireAttribute<double>("CrustThickness"));
    // Compute circumscribed radius
    m_circumscribedRadius = m_convex->computeCircumscribedRadius();
    // Volume and mass
    T volume  = m_convex->computeVolume();
    T density = T(0);
    m_mass    = T(0);
    if(element.hasAttribute("Density"))
    {
        density = T(element.requireAttribute<double>("Density"));
        m_mass  = density * volume;
    }
    // Equivalent radius: radius of sphere with same volume
    m_equivalentRadius = cbrt(T(3) * volume / (T(4) * PI<T>));
    // Material
    std::string material = element.requireAttribute<std::string>("Material");
    // checking if the material name is already defined. If yes, we access the ID and store it for
    // the rigid body. If it is not, we add the material to the map. Getting the ID of the last
    // material added to the map. This is basically the same as the size of the map.
    if(GrainsParameters<T>::m_materialMap.count(material) == 0)
    {
        uint id = GrainsParameters<T>::m_materialMap.size();
        GrainsParameters<T>::m_materialMap.emplace(material, id);
    }
    m_material = GrainsParameters<T>::m_materialMap[material];

    setInertia();
}
// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class RigidBody<float>;
template class RigidBody<double>;
