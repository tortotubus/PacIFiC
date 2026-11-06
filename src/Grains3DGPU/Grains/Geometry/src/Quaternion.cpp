#include "QuaternionImpl.hh"

// -------------------------------------------------------------------------------------------------
// Constructor from an XML node
template <typename T>
__HOST__ Quaternion<T>::Quaternion(grains::xml::OptionalElement root)
{
    const auto& element = grains::xml::requireElement(root);
    std::string type = element.requireAttribute<std::string>("Type");
    if(type == "Matrix")
    {
        Matrix3<T> mat(root);
        setQuaternion(mat);
    }
    else if(type == "Angles")
    {
        // read in radiands
        T aX = RADS_PER_DEG<T> * T(element.requireAttribute<double>("aX"));
        T aY = RADS_PER_DEG<T> * T(element.requireAttribute<double>("aY"));
        T aZ = RADS_PER_DEG<T> * T(element.requireAttribute<double>("aZ"));

        setQuaternion(aX, aY, aZ);
    }
    else if(type == "Identity")
    {
        setQuaternion(Vector3<T>(0, 0, 0), T(1));
    }
    else
        GAbort("A quaternion in one of the AngularPosition XML nodes is"
               " not a rotation matrix or angle.");
}
// -------------------------------------------------------------------------------------------------
// Output operator
template <typename T>
std::ostream& operator<<(std::ostream& fileOut, const Quaternion<T>& q)
{
    fileOut << q.getScalar() << " " << q.getVector();
    return (fileOut);
}
// -------------------------------------------------------------------------------------------------
// Input operator
template <typename T>
std::istream& operator>>(std::istream& fileIn, Quaternion<T>& q)
{
    Vector3<T> vec;
    T          scalar;
    fileIn >> scalar >> vec;
    q.setScalar(scalar);
    q.setVector(vec);
    return (fileIn);
}
// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class Quaternion<float>;
template class Quaternion<double>;

#define X(T)                                                                               \
    template std::ostream& operator<< <T>(std::ostream & fileOut, const Quaternion<T>& q); \
                                                                                           \
    template std::istream& operator>> <T>(std::istream & fileIn, Quaternion<T> & q);
X(float)
X(double)
#undef X
