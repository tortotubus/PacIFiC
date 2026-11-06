#include "ContactForceModelFactory.hh"
#include "GrainsParameters.hh"
#include "HookeContactForceModel.hh"
#include "HookeMemoryContactForceModel.hh"

/* ============================================================================================== */
/* High-Level Methods                                                                             */
/* ============================================================================================== */
// Creates and stores a ContactForceModel object in the host memory.
template <typename T>
__HOST__ void
    ContactForceModelFactory<T>::create(grains::xml::OptionalElement                           root,
                                        GrainsMemBuffer<ContactForceModel<T>*, MemType::HOST>& CF)
{
    uint numContactPairs = GrainsParameters<T>::m_numContactPairs;
    CF.initialize(numContactPairs);
    const auto allContacts
        = grains::xml::requireElement(root).getDescendants("ContactForceModel");
    for(const auto& contact : allContacts)
    {
        std::string contactType = contact.requireAttribute<std::string>("Type");
        // Contact Pair
        auto material = contact.getChild("Material");
        const auto& materialElement = grains::xml::requireElement(material);
        std::string matA = materialElement.requireAttribute<std::string>("materialA");
        std::string matB = materialElement.requireAttribute<std::string>("materialB");
        uint        matA_id  = GrainsParameters<T>::m_materialMap[matA];
        uint        matB_id  = GrainsParameters<T>::m_materialMap[matB];
        uint        index    = triangularHash(matA_id, matB_id);
        // Contact Parameters
        auto parameters = contact.getChild("Parameters");
        if(contactType == "Hooke")
            CF[index] = new HookeContactForceModel<T>(parameters);
        else if(contactType == "HookeMemory")
        {
            GrainsParameters<T>::m_isContactWithMemory = true;
            CF[index] = new HookeMemoryContactForceModel<T>(parameters);
        }
        else
            GAbort("Unknown contact force model! Aborting Grains!");
    }
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class ContactForceModelFactory<float>;
template class ContactForceModelFactory<double>;
