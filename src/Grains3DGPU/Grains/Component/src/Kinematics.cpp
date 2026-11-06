#include "KinematicsImpl.hh"

// -------------------------------------------------------------------------------------------------
// Output operator
template <typename T>
__HOST__ std::ostream& operator<<(std::ostream& fileOut, const Kinematics<T>& k)
{
    fileOut << "Translational Component: " << k.getTranslationalComponent() << "\n"
            << "Angular Component: " << k.getAngularComponent();
    return (fileOut);
}

// -------------------------------------------------------------------------------------------------
// Input operator
template <typename T>
__HOST__ std::istream& operator>>(std::istream& fileIn, Kinematics<T>& k)
{
    Vector3<T> vec;
    fileIn >> vec;
    k.setTranslationalComponent(vec);
    fileIn >> vec;
    k.setAngularComponent(vec);
    return (fileIn);
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class Kinematics<float>;
template class Kinematics<double>;

#define X(T)                                                                               \
    template std::ostream& operator<< <T>(std::ostream & fileOut, const Kinematics<T>& k); \
                                                                                           \
    template std::istream& operator>> <T>(std::istream & fileIn, Kinematics<T> & k);
X(float)
X(double)
#undef X
