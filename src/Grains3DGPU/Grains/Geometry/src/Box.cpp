#include "BoxImpl.hh"

// -------------------------------------------------------------------------------------------------
// Constructor with an input stream
template <typename T>
__HOST__ Box<T>::Box(std::istream& fileIn)
{
    readConvex(fileIn);
}
// -------------------------------------------------------------------------------------------------
// Constructor with an XML node as an input parameter
template <typename T>
__HOST__ Box<T>::Box(grains::xml::OptionalElement root)
{
    const auto& element = grains::xml::requireElement(root);
    m_extent[X] = T(element.requireAttribute<double>("LX")) / T(2);
    m_extent[Y] = T(element.requireAttribute<double>("LY")) / T(2);
    m_extent[Z] = T(element.requireAttribute<double>("LZ")) / T(2);
}
// -------------------------------------------------------------------------------------------------
// Input operator
template <typename T>
__HOST__ void Box<T>::readConvex(std::istream& fileIn)
{
    fileIn >> m_extent;
}
// -------------------------------------------------------------------------------------------------
// Output operator
template <typename T>
__HOST__ void Box<T>::writeConvex(std::ostream& fileOut) const
{
    fileOut << "Box: " << T(2) * m_extent << ".\n";
}
// -------------------------------------------------------------------------------------------------
// Returns the number of points to write the box in a Paraview format
template <typename T>
__HOST__ int Box<T>::numberOfPoints_PARAVIEW() const
{
    return (8);
}
// -------------------------------------------------------------------------------------------------
// Returns the number of elementary polytopes to write the box in a Paraview
// format
template <typename T>
__HOST__ int Box<T>::numberOfCells_PARAVIEW() const
{
    return (1);
}
// -------------------------------------------------------------------------------------------------
// Returns a list of points describing the box in a Paraview format
template <typename T>
__HOST__ std::list<Vector3<T>> Box<T>::writePoints_PARAVIEW(const Transform3<T>& transform,
                                                            Vector3<T> const*    translation) const
{
    std::list<Vector3<T>> ParaviewPoints;
    Vector3<T>            p;
    p.setValue(-m_extent[X], -m_extent[Y], -m_extent[Z]);
    ParaviewPoints.push_back(transform(p));
    p.setValue(m_extent[X], -m_extent[Y], -m_extent[Z]);
    ParaviewPoints.push_back(transform(p));
    p.setValue(m_extent[X], -m_extent[Y], m_extent[Z]);
    ParaviewPoints.push_back(transform(p));
    p.setValue(-m_extent[X], -m_extent[Y], m_extent[Z]);
    ParaviewPoints.push_back(transform(p));
    p.setValue(-m_extent[X], m_extent[Y], -m_extent[Z]);
    ParaviewPoints.push_back(transform(p));
    p.setValue(m_extent[X], m_extent[Y], -m_extent[Z]);
    ParaviewPoints.push_back(transform(p));
    p.setValue(m_extent[X], m_extent[Y], m_extent[Z]);
    ParaviewPoints.push_back(transform(p));
    p.setValue(-m_extent[X], m_extent[Y], m_extent[Z]);
    ParaviewPoints.push_back(transform(p));
    return (ParaviewPoints);
}
// -------------------------------------------------------------------------------------------------
// Writes the connectivity of the box in a Paraview format
template <typename T>
__HOST__ void Box<T>::writeConnection_PARAVIEW(std::list<uint>& connectivity,
                                               std::list<uint>& offsets,
                                               std::list<uint>& cellstype,
                                               uint&            firstpoint_globalnumber,
                                               uint&            last_offset) const
{
    uint count = firstpoint_globalnumber;
    for(uint i = 0; i < 8; ++i)
    {
        connectivity.push_back(count);
        ++count;
    }
    last_offset += 8;
    offsets.push_back(last_offset);
    cellstype.push_back(12);
    firstpoint_globalnumber += 8;
}
// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class Box<float>;
template class Box<double>;
