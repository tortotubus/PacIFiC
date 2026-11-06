#ifndef GRAINS_CONTACT_FORCE_MODEL_IMPL_HH
#define GRAINS_CONTACT_FORCE_MODEL_IMPL_HH

#include "ContactForceModel.hh"

// -------------------------------------------------------------------------------------------------
// Default constructor
template <typename T>
__HOSTDEVICE__ ContactForceModel<T>::ContactForceModel()
{
}

// -------------------------------------------------------------------------------------------------
// Copy constructor
template <typename T>
__HOSTDEVICE__ ContactForceModel<T>::ContactForceModel(ContactForceModel<T> const& cf)
{
}

// -------------------------------------------------------------------------------------------------
// Destructor
template <typename T>
__HOSTDEVICE__ ContactForceModel<T>::~ContactForceModel()
{
}

#endif
