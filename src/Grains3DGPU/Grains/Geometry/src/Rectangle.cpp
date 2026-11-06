#include "RectangleImpl.hh"

// -------------------------------------------------------------------------------------------------
// Constructor with an input stream
template <typename T>
__HOST__ Rectangle<T>::Rectangle(std::istream& fileIn)
{
    readConvex(fileIn);
}
// -------------------------------------------------------------------------------------------------
// Constructor with an XML node as an input parameter
template <typename T>
__HOST__ Rectangle<T>::Rectangle(grains::xml::OptionalElement root)
{
    const auto& element = grains::xml::requireElement(root);
    m_LX = T(element.requireAttribute<double>("LX")) / T(2);
    m_LY = T(element.requireAttribute<double>("LY")) / T(2);
}
// -------------------------------------------------------------------------------------------------
// Input operator
template <typename T>
__HOST__ void Rectangle<T>::readConvex(std::istream& fileIn)
{
    fileIn >> m_LX >> m_LY;
}
// -------------------------------------------------------------------------------------------------
// Output operator
template <typename T>
__HOST__ void Rectangle<T>::writeConvex(std::ostream& fileOut) const
{
    fileOut << "Rectangle: " << T(2) * m_LX << ", " << T(2) * m_LY << ".\n";
}
// -------------------------------------------------------------------------------------------------
// Returns the number of points to write the rectangle in a Paraview format
template <typename T>
__HOST__ int Rectangle<T>::numberOfPoints_PARAVIEW() const
{
    return (4);
}
// -------------------------------------------------------------------------------------------------
// Returns the number of elementary polytopes to write the rectangle in a
// Paraview format
template <typename T>
__HOST__ int Rectangle<T>::numberOfCells_PARAVIEW() const
{
    return (1);
}
// -------------------------------------------------------------------------------------------------
// Returns a list of points describing the rectangle in a Paraview format
template <typename T>
__HOST__ std::list<Vector3<T>>
         Rectangle<T>::writePoints_PARAVIEW(const Transform3<T>& transform,
                                       Vector3<T> const*    translation) const
{
    std::list<Vector3<T>> ParaviewPoints;
    Vector3<T>            p;
    p.setValue(-m_LX, -m_LY, T(0));
    ParaviewPoints.push_back(transform(p));
    p.setValue(m_LX, -m_LY, T(0));
    ParaviewPoints.push_back(transform(p));
    p.setValue(-m_LX, m_LY, T(0));
    ParaviewPoints.push_back(transform(p));
    p.setValue(m_LX, m_LY, T(0));
    ParaviewPoints.push_back(transform(p));
    return (ParaviewPoints);
}
// -------------------------------------------------------------------------------------------------
// Writes the connectivity of the box in a Paraview format
template <typename T>
__HOST__ void Rectangle<T>::writeConnection_PARAVIEW(std::list<uint>& connectivity,
                                                     std::list<uint>& offsets,
                                                     std::list<uint>& cellstype,
                                                     uint&            firstpoint_globalnumber,
                                                     uint&            last_offset) const
{
    uint count = firstpoint_globalnumber;
    for(uint i = 0; i < 4; ++i)
    {
        connectivity.push_back(count);
        ++count;
    }
    last_offset += 4;
    offsets.push_back(last_offset);
    cellstype.push_back(8);

    firstpoint_globalnumber += 4;
}
// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class Rectangle<float>;
template class Rectangle<double>;
