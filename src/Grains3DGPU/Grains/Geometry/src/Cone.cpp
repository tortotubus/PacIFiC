#include "ConeImpl.hh"

// -------------------------------------------------------------------------------------------------
// Constructor with an input stream
template <typename T>
__HOST__ Cone<T>::Cone(std::istream& fileIn)
{
    readConvex(fileIn);
}
// -------------------------------------------------------------------------------------------------
// Constructor with an XML node as an input parameter
template <typename T>
__HOST__ Cone<T>::Cone(grains::xml::OptionalElement root)
{
    const auto& element = grains::xml::requireElement(root);
    m_bottomRadius  = T(element.requireAttribute<double>("Radius"));
    m_quarterHeight = T(element.requireAttribute<double>("Height")) / T(4);
    m_sinAngle
        = m_bottomRadius
          / sqrt(m_bottomRadius * m_bottomRadius + T(16) * m_quarterHeight * m_quarterHeight);
}
// -------------------------------------------------------------------------------------------------
// Input operator
template <typename T>
__HOST__ void Cone<T>::readConvex(std::istream& fileIn)
{
    fileIn >> m_bottomRadius >> m_quarterHeight;
    m_quarterHeight /= T(4);
    m_sinAngle
        = m_bottomRadius
          / sqrt(m_bottomRadius * m_bottomRadius + T(16) * m_quarterHeight * m_quarterHeight);
}
// -------------------------------------------------------------------------------------------------
// Output operator
template <typename T>
__HOST__ void Cone<T>::writeConvex(std::ostream& fileOut) const
{
    fileOut << "Cone: " << m_bottomRadius << ", " << T(4) * m_quarterHeight << ".\n";
}
// -------------------------------------------------------------------------------------------------
// Returns the number of points to write the cone in a Paraview format
template <typename T>
__HOST__ int Cone<T>::numberOfPoints_PARAVIEW() const
{
    return (visuNodeNbOnPer + 2);
}
// -------------------------------------------------------------------------------------------------
// Returns the number of elementary polytopes to write the cone in a
// Paraview format
template <typename T>
__HOST__ int Cone<T>::numberOfCells_PARAVIEW() const
{
    return (visuNodeNbOnPer);
}
// -------------------------------------------------------------------------------------------------
// Returns a list of points describing the cone in a Paraview format
template <typename T>
__HOST__ std::list<Vector3<T>> Cone<T>::writePoints_PARAVIEW(const Transform3<T>& transform,
                                                             Vector3<T> const*    translation) const
{
    list<Vector3<T>> ParaviewPoints;
    Vector3<T>       pp, p;
    T                dtheta = TWO_PI<T> / visuNodeNbOnPer;

    // Disk rim
    p[Y] = -m_quarterHeight;
    for(int i = 0; i < visuNodeNbOnPer; ++i)
    {
        p[X] = m_bottomRadius * cos(i * dtheta);
        p[Z] = m_bottomRadius * sin(i * dtheta);
        pp   = transform(p);
        if(translation)
            pp += *translation;
        ParaviewPoints.push_back(pp);
    }

    // Disk center
    p[X] = T(0);
    p[Y] = -m_quarterHeight;
    p[Z] = T(0);
    pp   = transform(p);
    if(translation)
        pp += *translation;
    ParaviewPoints.push_back(pp);

    // Upper tip
    p[X] = T(0);
    p[Y] = T(3) * m_quarterHeight;
    p[Z] = T(0);
    pp   = transform(p);
    if(translation)
        pp += *translation;
    ParaviewPoints.push_back(pp);

    return (ParaviewPoints);
}
// -------------------------------------------------------------------------------------------------
// Writes the connectivity of the cone in a Paraview format
template <typename T>
__HOST__ void Cone<T>::writeConnection_PARAVIEW(std::list<uint>& connectivity,
                                                std::list<uint>& offsets,
                                                std::list<uint>& cellstype,
                                                uint&            firstpoint_globalnumber,
                                                uint&            last_offset) const
{
    for(uint i = 0; i < visuNodeNbOnPer - 1; ++i)
    {
        connectivity.push_back(firstpoint_globalnumber + i);
        connectivity.push_back(firstpoint_globalnumber + i + 1);
        connectivity.push_back(firstpoint_globalnumber + visuNodeNbOnPer);
        connectivity.push_back(firstpoint_globalnumber + visuNodeNbOnPer + 1);
        last_offset += 4;
        offsets.push_back(last_offset);
        cellstype.push_back(10);
    }
    connectivity.push_back(firstpoint_globalnumber + visuNodeNbOnPer - 1);
    connectivity.push_back(firstpoint_globalnumber);
    connectivity.push_back(firstpoint_globalnumber + visuNodeNbOnPer);
    connectivity.push_back(firstpoint_globalnumber + visuNodeNbOnPer + 1);
    last_offset += 4;
    offsets.push_back(last_offset);
    cellstype.push_back(10);

    firstpoint_globalnumber += visuNodeNbOnPer + 2;
}
// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class Cone<float>;
template class Cone<double>;

#undef visuNodeNbOnPer
