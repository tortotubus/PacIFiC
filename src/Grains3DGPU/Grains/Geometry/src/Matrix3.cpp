#include "Matrix3Impl.hh"

// -------------------------------------------------------------------------------------------------
// Constructor with an XML node
template <typename T>
__HOST__ Matrix3<T>::Matrix3(grains::xml::OptionalElement root)
{
    if(root)
    {
        std::string        values = root->text();
        std::istringstream inValues(values.c_str());
        inValues >> this->m_comp[XX] >> this->m_comp[XY] >> this->m_comp[XZ] >> this->m_comp[YX]
            >> this->m_comp[YY] >> this->m_comp[YZ] >> this->m_comp[ZX] >> this->m_comp[ZY]
            >> this->m_comp[ZZ];
    }
}
// -------------------------------------------------------------------------------------------------
// Output operator
template <typename T>
__HOST__ std::ostream& operator<<(std::ostream& fileOut, const Matrix3<T>& m)
{
    fileOut << m[X] << "\n" << m[Y] << "\n" << m[Z];
    return (fileOut);
}
// -------------------------------------------------------------------------------------------------
// Input operator
template <typename T>
__HOST__ std::istream& operator>>(std::istream& fileIn, Matrix3<T>& m)
{
    fileIn >> m[X] >> m[Y] >> m[Z];
    return (fileIn);
}
// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class Matrix3<float>;
template class Matrix3<double>;

#define X(T)                                                                            \
    template std::ostream& operator<< <T>(std::ostream & fileOut, const Matrix3<T>& m); \
                                                                                        \
    template std::istream& operator>> <T>(std::istream & fileIn, Matrix3<T> & m);
X(float)
X(double)
#undef X
