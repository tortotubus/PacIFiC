#include "ContactInfoImpl.hh"

// -------------------------------------------------------------------------------------------------
// Output operator
template <typename T>
__HOST__ std::ostream& operator<<(std::ostream& fileOut, const ContactInfo<T>& c)
{
    // Orientation first, followed by the position
    fileOut << "Contact Point: " << c.getContactPoint() << "\n"
            << "Contact Vector: " << c.getContactVector() << "\n"
            << "Overlap Distance: " << c.getOverlapDistance();
    return (fileOut);
}
// -------------------------------------------------------------------------------------------------
// Input operator
template <typename T>
__HOST__ std::istream& operator>>(std::istream& fileIn, ContactInfo<T>& c)
{
    GAbort("Input operator for ContactInfo is not implemented yet!");
    return (fileIn);
}
// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class ContactInfo<float>;
template class ContactInfo<double>;

#define X(T)                                                                                \
    template std::ostream& operator<< <T>(std::ostream & fileOut, const ContactInfo<T>& t); \
                                                                                            \
    template std::istream& operator>> <T>(std::istream & fileIn, ContactInfo<T> & t);
X(float)
X(double)
#undef X
