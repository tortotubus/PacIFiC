#include "CylinderImpl.hh"

// -------------------------------------------------------------------------------------------------
// Constructor with an input stream
template <typename T>
__HOST__ Cylinder<T>::Cylinder(std::istream& fileIn)
{
    readConvex(fileIn);
}
// -------------------------------------------------------------------------------------------------
// Constructor with an XML node as an input parameter
template <typename T>
__HOST__ Cylinder<T>::Cylinder(grains::xml::OptionalElement root)
{
    const auto& element = grains::xml::requireElement(root);
    m_radius     = T(element.requireAttribute<double>("Radius"));
    m_halfHeight = T(element.requireAttribute<double>("Height")) / T(2);
}
// -------------------------------------------------------------------------------------------------
// Input operator
template <typename T>
__HOST__ void Cylinder<T>::readConvex(std::istream& fileIn)
{
    fileIn >> m_radius >> m_halfHeight;
    m_halfHeight /= T(2);
}
// -------------------------------------------------------------------------------------------------
// Output operator
template <typename T>
__HOST__ void Cylinder<T>::writeConvex(std::ostream& fileOut) const
{
    fileOut << "Cylinder: " << m_radius << ", " << T(2) * m_halfHeight << ".\n";
}
// -------------------------------------------------------------------------------------------------
// Returns the number of points to write the cylinder in a Paraview format
template <typename T>
__HOST__ int Cylinder<T>::numberOfPoints_PARAVIEW() const
{
    return (2 * visuNodeNbOnPer + 2);
}
// -------------------------------------------------------------------------------------------------
// Returns the number of elementary polytopes to write the cylinder in a
// Paraview format
template <typename T>
__HOST__ int Cylinder<T>::numberOfCells_PARAVIEW() const
{
    return (visuNodeNbOnPer);
}
// -------------------------------------------------------------------------------------------------
// Returns a list of points describing the cylinder in a Paraview format
template <typename T>
__HOST__ std::list<Vector3<T>>
         Cylinder<T>::writePoints_PARAVIEW(const Transform3<T>& transform,
                                      Vector3<T> const*    translation) const
{
    list<Vector3<T>> ParaviewPoints;
    Vector3<T>       pp, p;
    T                dtheta = TWO_PI<T> / visuNodeNbOnPer;

    // Lower disk rim
    p[Y] = -m_halfHeight;
    for(int i = 0; i < visuNodeNbOnPer; ++i)
    {
        p[X] = m_radius * cos(i * dtheta);
        p[Z] = m_radius * sin(i * dtheta);
        pp   = transform(p);
        if(translation)
            pp += *translation;
        ParaviewPoints.push_back(pp);
    }

    // Upper disk rim
    p[Y] = m_halfHeight;
    for(int i = 0; i < visuNodeNbOnPer; ++i)
    {
        p[X] = m_radius * cos(i * dtheta);
        p[Z] = m_radius * sin(i * dtheta);
        pp   = transform(p);
        if(translation)
            pp += *translation;
        ParaviewPoints.push_back(pp);
    }

    // Lower disk center
    p[X] = T(0);
    p[Y] = -m_halfHeight;
    p[Z] = T(0);
    pp   = transform(p);
    if(translation)
        pp += *translation;
    ParaviewPoints.push_back(pp);

    // Upper disk center
    p[Y] = m_halfHeight;
    pp   = transform(p);
    if(translation)
        pp += *translation;
    ParaviewPoints.push_back(pp);

    return (ParaviewPoints);
}
// -------------------------------------------------------------------------------------------------
// Writes the connectivity of the cylinder in a Paraview format
template <typename T>
__HOST__ void Cylinder<T>::writeConnection_PARAVIEW(std::list<uint>& connectivity,
                                                    std::list<uint>& offsets,
                                                    std::list<uint>& cellstype,
                                                    uint&            firstpoint_globalnumber,
                                                    uint&            last_offset) const
{
    for(int i = 0; i < visuNodeNbOnPer - 1; ++i)
    {
        connectivity.push_back(firstpoint_globalnumber + i);
        connectivity.push_back(firstpoint_globalnumber + i + 1);
        connectivity.push_back(firstpoint_globalnumber + 2 * visuNodeNbOnPer);
        connectivity.push_back(firstpoint_globalnumber + i + visuNodeNbOnPer);
        connectivity.push_back(firstpoint_globalnumber + i + visuNodeNbOnPer + 1);
        connectivity.push_back(firstpoint_globalnumber + 2 * visuNodeNbOnPer + 1);
        last_offset += 6;
        offsets.push_back(last_offset);
        cellstype.push_back(13);
    }
    connectivity.push_back(firstpoint_globalnumber + visuNodeNbOnPer - 1);
    connectivity.push_back(firstpoint_globalnumber);
    connectivity.push_back(firstpoint_globalnumber + 2 * visuNodeNbOnPer);
    connectivity.push_back(firstpoint_globalnumber + 2 * visuNodeNbOnPer - 1);
    connectivity.push_back(firstpoint_globalnumber + visuNodeNbOnPer);
    connectivity.push_back(firstpoint_globalnumber + 2 * visuNodeNbOnPer + 1);
    last_offset += 6;
    offsets.push_back(last_offset);
    cellstype.push_back(13);

    firstpoint_globalnumber += 2 * visuNodeNbOnPer + 2;
}
// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class Cylinder<float>;
template class Cylinder<double>;

#undef visuNodeNbOnPer
