#include "GrainsPlatform.hh"
#include "Vector3.hh"
#include "VectorMath.hh"

// -------------------------------------------------------------------------------------------------
// Constructor from an XML node
template <typename T>
__HOST__ Vector3<T>::Vector3(grains::xml::OptionalElement root)
{
    const auto& element = grains::xml::requireElement(root);
    T x = T(element.requireAttribute<double>("X"));
    T y = T(element.requireAttribute<double>("Y"));
    T z = T(element.requireAttribute<double>("Z"));
    setValue(x, y, z);
}

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
// Output operator
template <typename T>
__HOST__ std::ostream& operator<<(std::ostream& fileOut, const Vector3<T>& v)
{
    fileOut << "[" << v[X] << ", " << v[Y] << ", " << v[Z] << "]";
    return (fileOut);
}

// -------------------------------------------------------------------------------------------------
// Input operator
template <typename T>
__HOST__ std::istream& operator>>(std::istream& fileIn, Vector3<T>& v)
{
    fileIn >> v[X] >> v[Y] >> v[Z];
    return (fileIn);
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class Vector3<float>;
template class Vector3<double>;

#define X(T)                                                                            \
    template std::ostream& operator<< <T>(std::ostream & fileOut, const Vector3<T>& v); \
                                                                                        \
    template std::istream& operator>> <T>(std::istream & fileIn, Vector3<T> & v);
X(float)
X(double)
#undef X
