#include "Octahedron.hh"
#include "VectorMath.hh"

namespace
{
    template <typename T>
    __HOSTDEVICE__ Vector3<T> octahedronVertex(int idx, T radius)
    {
        switch(idx)
        {
        case 0:
            return Vector3<T>(radius, T(0), T(0));
        case 1:
            return Vector3<T>(-radius, T(0), T(0));
        case 2:
            return Vector3<T>(T(0), radius, T(0));
        case 3:
            return Vector3<T>(T(0), -radius, T(0));
        case 4:
            return Vector3<T>(T(0), T(0), radius);
        default:
            return Vector3<T>(T(0), T(0), -radius);
        }
    }

    constexpr int kOctahedronFaces[8][3]
        = {{0, 4, 2}, {2, 4, 1}, {1, 4, 3}, {3, 4, 0}, {0, 2, 5}, {2, 1, 5}, {1, 3, 5}, {3, 0, 5}};
}

// -------------------------------------------------------------------------------------------------
// Constructor with circumradius as input parameter
template <typename T>
__HOSTDEVICE__ Octahedron<T>::Octahedron(T radius)
    : m_radius(radius)
{
}

// -------------------------------------------------------------------------------------------------
// Constructor with an input stream
template <typename T>
__HOST__ Octahedron<T>::Octahedron(std::istream& fileIn)
{
    readConvex(fileIn);
}

// -------------------------------------------------------------------------------------------------
// Constructor with an XML node as an input parameter
template <typename T>
__HOST__ Octahedron<T>::Octahedron(DOMNode* root)
{
    m_radius = T(ReaderXML::getNodeAttr_Double(root, "Radius"));
}

// -------------------------------------------------------------------------------------------------
// Destructor
template <typename T>
__HOSTDEVICE__ Octahedron<T>::~Octahedron()
{
}

// -------------------------------------------------------------------------------------------------
// Returns the convex type
template <typename T>
__HOSTDEVICE__ ConvexType Octahedron<T>::getConvexType() const
{
    return OCTAHEDRON;
}

// -------------------------------------------------------------------------------------------------
// Gets the circumradius
template <typename T>
__HOSTDEVICE__ T Octahedron<T>::getRadius() const
{
    return m_radius;
}

// -------------------------------------------------------------------------------------------------
// Returns a clone of the octahedron
template <typename T>
__HOSTDEVICE__ Convex<T>* Octahedron<T>::clone() const
{
    return new Octahedron<T>(m_radius);
}

// -------------------------------------------------------------------------------------------------
// Returns the volume of the octahedron
template <typename T>
__HOSTDEVICE__ T Octahedron<T>::computeVolume() const
{
    return T(4) * m_radius * m_radius * m_radius / T(3);
}

// -------------------------------------------------------------------------------------------------
// Computes the diagonal inertia tensor
template <typename T>
__HOSTDEVICE__ void Octahedron<T>::computeInertia(T (&inertia)[3]) const
{
    const T geometricInertia = computeVolume() * m_radius * m_radius / T(5);
    inertia[0]               = geometricInertia;
    inertia[1]               = geometricInertia;
    inertia[2]               = geometricInertia;
}

// -------------------------------------------------------------------------------------------------
// Returns the circumscribed radius of the octahedron
template <typename T>
__HOSTDEVICE__ T Octahedron<T>::computeCircumscribedRadius() const
{
    return m_radius;
}

// -------------------------------------------------------------------------------------------------
// Returns the bounding box to octahedron
template <typename T>
__HOSTDEVICE__ Vector3<T> Octahedron<T>::computeBoundingBox() const
{
    return Vector3<T>(m_radius, m_radius, m_radius);
}

// -------------------------------------------------------------------------------------------------
// Returns the bounding cylinder to octahedron
template <typename T>
__HOSTDEVICE__ Vector3<T> Octahedron<T>::computeBoundingCylinder() const
{
    return Vector3<T>(m_radius, m_radius, T(1));
}

// -------------------------------------------------------------------------------------------------
// Octahedron support function, returns the support point P, i.e. the point on the
// surface of the octahedron that satisfies max(P.v)
template <typename T>
__HOSTDEVICE__ Vector3<T> Octahedron<T>::support(const Vector3<T>& v) const
{
    const T ax = fabs(v[X]);
    const T ay = fabs(v[Y]);
    const T az = fabs(v[Z]);
    if(ax >= ay && ax >= az)
        return Vector3<T>(v[X] < T(0) ? -m_radius : m_radius, T(0), T(0));
    if(ay >= az)
        return Vector3<T>(T(0), v[Y] < T(0) ? -m_radius : m_radius, T(0));
    return Vector3<T>(T(0), T(0), v[Z] < T(0) ? -m_radius : m_radius);
}

// -------------------------------------------------------------------------------------------------
// Eroded octahedron support: evaluates support of octahedron with radius shrunk by crust
template <typename T>
__HOSTDEVICE__ Vector3<T> Octahedron<T>::support(const Vector3<T>& v, T crust) const
{
    const T radius = m_radius - crust;
    const T ax     = fabs(v[X]);
    const T ay     = fabs(v[Y]);
    const T az     = fabs(v[Z]);
    if(ax >= ay && ax >= az)
        return Vector3<T>(v[X] < T(0) ? -radius : radius, T(0), T(0));
    if(ay >= az)
        return Vector3<T>(T(0), v[Y] < T(0) ? -radius : radius, T(0));
    return Vector3<T>(T(0), T(0), v[Z] < T(0) ? -radius : radius);
}

// -------------------------------------------------------------------------------------------------
// Input operator
template <typename T>
__HOST__ void Octahedron<T>::readConvex(std::istream& fileIn)
{
    fileIn >> m_radius;
}

// -------------------------------------------------------------------------------------------------
// Output operator
template <typename T>
__HOST__ void Octahedron<T>::writeConvex(std::ostream& fileOut) const
{
    fileOut << "Octahedron: " << m_radius << ".\n";
}

// -------------------------------------------------------------------------------------------------
// Returns the number of points to write the octahedron in a Paraview format
template <typename T>
__HOST__ int Octahedron<T>::numberOfPoints_PARAVIEW() const
{
    return 6;
}

// -------------------------------------------------------------------------------------------------
// Returns the number of elementary polytopes to write the octahedron in a Paraview format
template <typename T>
__HOST__ int Octahedron<T>::numberOfCells_PARAVIEW() const
{
    return 8;
}

// -------------------------------------------------------------------------------------------------
// Returns a list of points describing the octahedron in a Paraview format
template <typename T>
__HOST__ std::list<Vector3<T>>
         Octahedron<T>::writePoints_PARAVIEW(const Transform3<T>& transform,
                                        Vector3<T> const*    translation) const
{
    std::list<Vector3<T>> points;
    for(int i = 0; i < 6; ++i)
    {
        Vector3<T> out = transform(octahedronVertex<T>(i, m_radius));
        if(translation)
            out += *translation;
        points.push_back(out);
    }
    return points;
}

// -------------------------------------------------------------------------------------------------
// Writes the connectivity of the octahedron in a Paraview format
template <typename T>
__HOST__ void Octahedron<T>::writeConnection_PARAVIEW(std::list<uint>& connectivity,
                                                      std::list<uint>& offsets,
                                                      std::list<uint>& cellstype,
                                                      uint&            firstpoint_globalnumber,
                                                      uint&            last_offset) const
{
    for(const auto& face : kOctahedronFaces)
    {
        connectivity.push_back(firstpoint_globalnumber + static_cast<uint>(face[0]));
        connectivity.push_back(firstpoint_globalnumber + static_cast<uint>(face[1]));
        connectivity.push_back(firstpoint_globalnumber + static_cast<uint>(face[2]));
        last_offset += 3;
        offsets.push_back(last_offset);
        cellstype.push_back(5);
    }
    firstpoint_globalnumber += 6;
}

// -------------------------------------------------------------------------------------------------
// Returns the shape parameters
template <typename T>
__HOSTDEVICE__ void Octahedron<T>::getShapeParameters(T (&params)[5]) const
{
    params[0] = m_radius;
    params[1] = T(0);
    params[2] = T(0);
    params[3] = T(0);
    params[4] = T(0);
}

// -------------------------------------------------------------------------------------------------
// Returns the convex name
template <typename T>
__HOST__ std::string Octahedron<T>::getConvexName() const
{
    return "Octahedron";
}

// -------------------------------------------------------------------------------------------------
// Writes the octahedron in OBJ format
template <typename T>
__HOST__ void Octahedron<T>::writeOBJ(std::ostream& f, size_t& firstpoint_number) const
{
    for(int i = 0; i < 6; ++i)
    {
        const Vector3<T> vertex = octahedronVertex<T>(i, m_radius);
        f << "v " << vertex[X] << " " << vertex[Y] << " " << vertex[Z] << "\n";
    }

    for(const auto& face : kOctahedronFaces)
    {
        f << "f " << firstpoint_number + static_cast<size_t>(face[0]) << " "
          << firstpoint_number + static_cast<size_t>(face[1]) << " "
          << firstpoint_number + static_cast<size_t>(face[2]) << "\n";
    }

    firstpoint_number += 6;
}

// -------------------------------------------------------------------------------------------------
// Explicit template instantiations
template class Octahedron<float>;
template class Octahedron<double>;