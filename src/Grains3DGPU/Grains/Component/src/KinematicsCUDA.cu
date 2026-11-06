#include "KinematicsImpl.hh"

// Emit the device-callable float and double specializations. The ordinary C++ translation unit
// emits the corresponding host specializations and owns the host-only stream operators.
template class Kinematics<float>;
template class Kinematics<double>;
