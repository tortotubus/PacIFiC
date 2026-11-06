#include "GrainsPlatform.hh"
#include "Torce.hh"
#include "Vector3.hh"
#include "VectorMath.hh"

// -------------------------------------------------------------------------------------------------
// Output operator
template <typename T>
__HOST__ std::ostream& operator<<(std::ostream& fileOut, const Torce<T>& t)
{
    fileOut << "Torque: " << t.getTorque() << "\n"
            << "Force: " << t.getForce();
    return (fileOut);
}

// -------------------------------------------------------------------------------------------------
// Input operator
template <typename T>
__HOST__ std::istream& operator>>(std::istream& fileIn, Torce<T>& t)
{
    Vector3<T> vec;
    fileIn >> vec;
    t.setTorque(vec);
    fileIn >> vec;
    t.setForce(vec);
    return (fileIn);
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class Torce<float>;
template class Torce<double>;

#define X(T)                                                                          \
    template std::ostream& operator<< <T>(std::ostream & fileOut, const Torce<T>& t); \
    template std::istream& operator>> <T>(std::istream & fileIn, Torce<T> & t);
X(float)
X(double)
#undef X
