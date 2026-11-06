#include "SuperquadricImpl.hh"

// -------------------------------------------------------------------------------------------------
// Constructor with an input stream
template <typename T>
__HOST__ Superquadric<T>::Superquadric(std::istream& fileIn)
{
    readConvex(fileIn);
}
// -------------------------------------------------------------------------------------------------
// Constructor with an XML node as an input parameter
template <typename T>
__HOST__ Superquadric<T>::Superquadric(grains::xml::OptionalElement root)
{
    const auto& element = grains::xml::requireElement(root);
    m_a  = T(element.requireAttribute<double>("a"));
    m_b  = T(element.requireAttribute<double>("b"));
    m_c  = T(element.requireAttribute<double>("c"));
    m_n1 = T(element.requireAttribute<double>("n1"));
    m_n2 = T(element.requireAttribute<double>("n2"));
}
// -------------------------------------------------------------------------------------------------
// Input operator
template <typename T>
__HOST__ void Superquadric<T>::readConvex(std::istream& fileIn)
{
    std::cout << "Program Error :\n"
              << "Superquadric::readConvex is not implemented.\n";
    exit(3);
}
// -------------------------------------------------------------------------------------------------
// Output operator
template <typename T>
__HOST__ void Superquadric<T>::writeConvex(std::ostream& fileOut) const
{
    fileOut << "Superquadric: " << m_a << ", " << m_b << ", " << m_c << ", " << m_n1 << ", " << m_n2
            << ".\n";
}
// -------------------------------------------------------------------------------------------------
// Returns the number of points to write the superquadric in a Paraview format
template <typename T>
__HOST__ int Superquadric<T>::numberOfPoints_PARAVIEW() const
{
    return (visuNodeNbOnPer * (visuNodeNbOnPer - 1) + 3);
}
// -------------------------------------------------------------------------------------------------
// Returns the number of elementary polytopes to write the superquadric in a
// Paraview format
template <typename T>
__HOST__ int Superquadric<T>::numberOfCells_PARAVIEW() const
{
    return (visuNodeNbOnPer * visuNodeNbOnPer);
}
// -------------------------------------------------------------------------------------------------
// Returns a list of points describing the superquadric in a Paraview format
template <typename T>
__HOST__ std::list<Vector3<T>>
         Superquadric<T>::writePoints_PARAVIEW(const Transform3<T>& transform,
                                          Vector3<T> const*    translation) const
{
    T const               eps1 = T(2) / m_n1;
    T const               eps2 = T(2) / m_n2;
    std::list<Vector3<T>> ParaviewPoints;
    T const               dtheta = PI<T> / visuNodeNbOnPer;
    T const               dphi   = TWO_PI<T> / visuNodeNbOnPer;
    Vector3<T>            pp, pptrans;

    pp[X] = pp[Y] = T(0);
    // Gravity center
    pp[Z]   = T(0);
    pptrans = transform(pp);
    if(translation)
        pptrans += *translation;
    ParaviewPoints.push_back(pptrans);

    // Top point
    pp[Z]   = m_c;
    pptrans = transform(pp);
    if(translation)
        pptrans += *translation;
    ParaviewPoints.push_back(pptrans);

    // Bottom point
    pp[Z]   = -m_c;
    pptrans = transform(pp);
    if(translation)
        pptrans += *translation;
    ParaviewPoints.push_back(pptrans);

    // Regular points on the surface
    T cost, sint, costeps1, sinteps1, cosp, sinp;

    for(int i = 1; i < visuNodeNbOnPer; i++)
    {
        cost = cos(i * dtheta);
        sint = sin(i * dtheta);

        if(cost == T(0))
            costeps1 = T(0);
        else if(cost < T(0))
            costeps1 = -pow(-cost, eps1);
        else
            costeps1 = pow(cost, eps1);

        // Theta is always strictly between 0 and pi so sint is strictly
        // positive
        sinteps1 = pow(sint, eps1);

        for(int j = 0; j < visuNodeNbOnPer; j++)
        {
            cosp = cos(j * dphi);
            sinp = sin(j * dphi);
            if(cosp == T(0))
                pp[X] = T(0);
            else if(cosp < T(0))
                pp[X] = -m_a * sinteps1 * pow(-cosp, eps2);
            else
                pp[X] = m_a * sinteps1 * pow(cosp, eps2);

            if(sinp == T(0))
                pp[Y] = T(0);
            else if(sinp < T(0))
                pp[Y] = -m_b * sinteps1 * pow(-sinp, eps2);
            else
                pp[Y] = m_b * sinteps1 * pow(sinp, eps2);

            pp[Z] = m_c * costeps1;

            pptrans = transform(pp);
            if(translation)
                pptrans += *translation;
            ParaviewPoints.push_back(pptrans);
        }
    }

    return (ParaviewPoints);
}
// -------------------------------------------------------------------------------------------------
// Writes the connectivity of the superquadric in a Paraview format
template <typename T>
__HOST__ void Superquadric<T>::writeConnection_PARAVIEW(std::list<uint>& connectivity,
                                                        std::list<uint>& offsets,
                                                        std::list<uint>& cellstype,
                                                        uint&            firstpoint_globalnumber,
                                                        uint&            last_offset) const
{
    // Top cells: tetrahedron
    for(uint j = 0; j < visuNodeNbOnPer - 1; ++j)
    {
        // Center
        connectivity.push_back(firstpoint_globalnumber);
        // Top point
        connectivity.push_back(firstpoint_globalnumber + 1);
        connectivity.push_back(firstpoint_globalnumber + 3 + j);
        connectivity.push_back(firstpoint_globalnumber + 3 + j + 1);
        last_offset += 4;
        offsets.push_back(last_offset);
        cellstype.push_back(10);
    }
    // Center
    connectivity.push_back(firstpoint_globalnumber);
    // Top point
    connectivity.push_back(firstpoint_globalnumber + 1);
    // Last point of the first latitude (i=1)
    connectivity.push_back(firstpoint_globalnumber + 3 + visuNodeNbOnPer - 1);
    // First point of the first latitude (i=1)
    connectivity.push_back(firstpoint_globalnumber + 3);
    last_offset += 4;
    offsets.push_back(last_offset);
    cellstype.push_back(10);

    // Regular cells: Pyramid
    uint pointact = 3;  // Current point (we will browse the grid)

    for(uint i = 1; i < visuNodeNbOnPer - 1; ++i)
    {
        for(uint j = 0; j < visuNodeNbOnPer - 1; ++j)
        {
            // Center
            connectivity.push_back(firstpoint_globalnumber);
            // Current point
            connectivity.push_back(firstpoint_globalnumber + pointact);
            // same latitude, +1 longitude
            connectivity.push_back(firstpoint_globalnumber + pointact + 1);
            // +1 latitude, +1 longitude
            connectivity.push_back(firstpoint_globalnumber + pointact + visuNodeNbOnPer + 1);
            // +1 latitude, same longitude
            connectivity.push_back(firstpoint_globalnumber + pointact + visuNodeNbOnPer);
            last_offset += 5;
            offsets.push_back(last_offset);
            cellstype.push_back(14);
            pointact++;
        }
        // Center
        connectivity.push_back(firstpoint_globalnumber);
        // Current point (last of its latitude)
        connectivity.push_back(firstpoint_globalnumber + pointact);
        // First point (j=0) on the same latitude as pointact
        connectivity.push_back(firstpoint_globalnumber + pointact - (visuNodeNbOnPer - 1));
        // First point (j=0) on the latitude under that of pointact
        connectivity.push_back(firstpoint_globalnumber + pointact + 1);
        // +1 latitude, same longitude
        connectivity.push_back(firstpoint_globalnumber + pointact + visuNodeNbOnPer);
        last_offset += 5;
        offsets.push_back(last_offset);
        cellstype.push_back(14);
        pointact++;
    }

    // Bottom cells: tetrahedron
    const uint Firstptlastlat = 3 + visuNodeNbOnPer * (visuNodeNbOnPer - 2);

    for(uint j = 0; j < visuNodeNbOnPer - 1; ++j)
    {
        connectivity.push_back(firstpoint_globalnumber);      // Center
        connectivity.push_back(firstpoint_globalnumber + 2);  // Bottom point
        connectivity.push_back(firstpoint_globalnumber + Firstptlastlat + j);
        connectivity.push_back(firstpoint_globalnumber + Firstptlastlat + j + 1);
        last_offset += 4;
        offsets.push_back(last_offset);
        cellstype.push_back(10);
    }
    connectivity.push_back(firstpoint_globalnumber);      // Center
    connectivity.push_back(firstpoint_globalnumber + 2);  // Bottom point
    // Last point last latitude
    connectivity.push_back(firstpoint_globalnumber + Firstptlastlat + visuNodeNbOnPer - 1);
    // First point last latitude
    connectivity.push_back(firstpoint_globalnumber + Firstptlastlat);
    last_offset += 4;
    offsets.push_back(last_offset);
    cellstype.push_back(10);

    firstpoint_globalnumber += visuNodeNbOnPer * (visuNodeNbOnPer - 1) + 3;
}
// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class Superquadric<float>;
template class Superquadric<double>;
