#include "Dodecahedron.hh"
#include "VectorMath.hh"

namespace
{
    constexpr int kDodecahedronFaces[12][5] = {{0, 12, 1, 17, 16},
                                               {3, 13, 2, 16, 17},
                                               {4, 14, 12, 0, 8},
                                               {5, 19, 7, 11, 9},
                                               {6, 18, 4, 8, 10},
                                               {7, 15, 13, 3, 11},
                                               {9, 11, 3, 17, 1},
                                               {10, 8, 0, 16, 2},
                                               {14, 5, 9, 1, 12},
                                               {15, 6, 10, 2, 13},
                                               {18, 19, 5, 14, 4},
                                               {19, 18, 6, 15, 7}};

    template <typename T>
    __HOSTDEVICE__ T goldenRatio()
    {
        return (T(1) + sqrt(T(5))) / T(2);
    }

    template <typename T>
    __HOSTDEVICE__ T inverseGoldenRatio()
    {
        return (sqrt(T(5)) - T(1)) / T(2);
    }

    template <typename T>
    __HOSTDEVICE__ Vector3<T> dodecahedronVertex(int idx, T radius)
    {
        const T phi    = goldenRatio<T>();
        const T invPhi = inverseGoldenRatio<T>();
        const T scale  = radius / sqrt(T(3));
        switch(idx)
        {
        case 0:
            return Vector3<T>(-scale, -scale, -scale);
        case 1:
            return Vector3<T>(-scale, -scale, scale);
        case 2:
            return Vector3<T>(-scale, scale, -scale);
        case 3:
            return Vector3<T>(-scale, scale, scale);
        case 4:
            return Vector3<T>(scale, -scale, -scale);
        case 5:
            return Vector3<T>(scale, -scale, scale);
        case 6:
            return Vector3<T>(scale, scale, -scale);
        case 7:
            return Vector3<T>(scale, scale, scale);
        case 8:
            return Vector3<T>(T(0), -invPhi * scale, -phi * scale);
        case 9:
            return Vector3<T>(T(0), -invPhi * scale, phi * scale);
        case 10:
            return Vector3<T>(T(0), invPhi * scale, -phi * scale);
        case 11:
            return Vector3<T>(T(0), invPhi * scale, phi * scale);
        case 12:
            return Vector3<T>(-invPhi * scale, -phi * scale, T(0));
        case 13:
            return Vector3<T>(-invPhi * scale, phi * scale, T(0));
        case 14:
            return Vector3<T>(invPhi * scale, -phi * scale, T(0));
        case 15:
            return Vector3<T>(invPhi * scale, phi * scale, T(0));
        case 16:
            return Vector3<T>(-phi * scale, T(0), -invPhi * scale);
        case 17:
            return Vector3<T>(-phi * scale, T(0), invPhi * scale);
        case 18:
            return Vector3<T>(phi * scale, T(0), -invPhi * scale);
        default:
            return Vector3<T>(phi * scale, T(0), invPhi * scale);
        }
    }
}

// -------------------------------------------------------------------------------------------------
// Constructor with circumradius as input parameter
template <typename T>
__HOSTDEVICE__ Dodecahedron<T>::Dodecahedron(T radius)
    : m_radius(radius)
{
}

// -------------------------------------------------------------------------------------------------
// Constructor with an input stream
template <typename T>
__HOST__ Dodecahedron<T>::Dodecahedron(std::istream& fileIn)
{
    readConvex(fileIn);
}

// -------------------------------------------------------------------------------------------------
// Constructor with an XML node as an input parameter
template <typename T>
__HOST__ Dodecahedron<T>::Dodecahedron(DOMNode* root)
{
    m_radius = T(ReaderXML::getNodeAttr_Double(root, "Radius"));
}

// -------------------------------------------------------------------------------------------------
// Destructor
template <typename T>
__HOSTDEVICE__ Dodecahedron<T>::~Dodecahedron()
{
}

// -------------------------------------------------------------------------------------------------
// Returns the convex type
template <typename T>
__HOSTDEVICE__ ConvexType Dodecahedron<T>::getConvexType() const
{
    return DODECAHEDRON;
}

// -------------------------------------------------------------------------------------------------
// Gets the circumradius
template <typename T>
__HOSTDEVICE__ T Dodecahedron<T>::getRadius() const
{
    return m_radius;
}

// -------------------------------------------------------------------------------------------------
// Returns a clone of the dodecahedron
template <typename T>
__HOSTDEVICE__ Convex<T>* Dodecahedron<T>::clone() const
{
    return new Dodecahedron<T>(m_radius);
}

// -------------------------------------------------------------------------------------------------
// Returns the volume of the dodecahedron
template <typename T>
__HOSTDEVICE__ T Dodecahedron<T>::computeVolume() const
{
    return T(2) * (T(5) + sqrt(T(5))) * m_radius * m_radius * m_radius / (T(3) * sqrt(T(3)));
}

// -------------------------------------------------------------------------------------------------
// Computes the diagonal inertia tensor
template <typename T>
__HOSTDEVICE__ void Dodecahedron<T>::computeInertia(T (&inertia)[3]) const
{
    const T coeff            = (T(45) + T(11) * sqrt(T(5))) / T(225);
    const T geometricInertia = computeVolume() * coeff * m_radius * m_radius;
    inertia[0]               = geometricInertia;
    inertia[1]               = geometricInertia;
    inertia[2]               = geometricInertia;
}

// -------------------------------------------------------------------------------------------------
// Returns the circumscribed radius of the dodecahedron
template <typename T>
__HOSTDEVICE__ T Dodecahedron<T>::computeCircumscribedRadius() const
{
    return m_radius;
}

// -------------------------------------------------------------------------------------------------
// Returns the bounding box to dodecahedron
template <typename T>
__HOSTDEVICE__ Vector3<T> Dodecahedron<T>::computeBoundingBox() const
{
    const T extent = goldenRatio<T>() * m_radius / sqrt(T(3));
    return Vector3<T>(extent, extent, extent);
}

// -------------------------------------------------------------------------------------------------
// Returns the bounding cylinder to dodecahedron
template <typename T>
__HOSTDEVICE__ Vector3<T> Dodecahedron<T>::computeBoundingCylinder() const
{
    return Vector3<T>(m_radius, goldenRatio<T>() * m_radius / sqrt(T(3)), T(1));
}

// -------------------------------------------------------------------------------------------------
// Dodecahedron support function, returns the support point P, i.e. the point on the
// surface of the dodecahedron that satisfies max(P.v)
template <typename T>
__HOSTDEVICE__ Vector3<T> Dodecahedron<T>::support(const Vector3<T>& v) const
{
    Vector3<T> best      = dodecahedronVertex<T>(0, m_radius);
    T          bestScore = best * v;
    for(int i = 1; i < 20; ++i)
    {
        const Vector3<T> candidate = dodecahedronVertex<T>(i, m_radius);
        const T          score     = candidate * v;
        if(score > bestScore)
        {
            best      = candidate;
            bestScore = score;
        }
    }
    return best;
}

// -------------------------------------------------------------------------------------------------
// Eroded dodecahedron support: evaluates support of dodecahedron with radius shrunk by crust
template <typename T>
__HOSTDEVICE__ Vector3<T> Dodecahedron<T>::support(const Vector3<T>& v, T crust) const
{
    Vector3<T> best      = dodecahedronVertex<T>(0, m_radius - crust);
    T          bestScore = best * v;
    for(int i = 1; i < 20; ++i)
    {
        const Vector3<T> candidate = dodecahedronVertex<T>(i, m_radius - crust);
        const T          score     = candidate * v;
        if(score > bestScore)
        {
            best      = candidate;
            bestScore = score;
        }
    }
    return best;
}

// -------------------------------------------------------------------------------------------------
// Input operator
template <typename T>
__HOST__ void Dodecahedron<T>::readConvex(std::istream& fileIn)
{
    fileIn >> m_radius;
}

// -------------------------------------------------------------------------------------------------
// Output operator
template <typename T>
__HOST__ void Dodecahedron<T>::writeConvex(std::ostream& fileOut) const
{
    fileOut << "Dodecahedron: " << m_radius << ".\n";
}

// -------------------------------------------------------------------------------------------------
// Returns the number of points to write the dodecahedron in a Paraview format
template <typename T>
__HOST__ int Dodecahedron<T>::numberOfPoints_PARAVIEW() const
{
    return 20;
}

// -------------------------------------------------------------------------------------------------
// Returns the number of elementary polytopes to write the dodecahedron in a Paraview format
template <typename T>
__HOST__ int Dodecahedron<T>::numberOfCells_PARAVIEW() const
{
    return 12;
}

// -------------------------------------------------------------------------------------------------
// Returns a list of points describing the dodecahedron in a Paraview format
template <typename T>
__HOST__ std::list<Vector3<T>>
         Dodecahedron<T>::writePoints_PARAVIEW(const Transform3<T>& transform,
                                          Vector3<T> const*    translation) const
{
    std::list<Vector3<T>> points;
    for(int i = 0; i < 20; ++i)
    {
        Vector3<T> out = transform(dodecahedronVertex<T>(i, m_radius));
        if(translation)
            out += *translation;
        points.push_back(out);
    }
    return points;
}

// -------------------------------------------------------------------------------------------------
// Writes the connectivity of the dodecahedron in a Paraview format
template <typename T>
__HOST__ void Dodecahedron<T>::writeConnection_PARAVIEW(std::list<uint>& connectivity,
                                                        std::list<uint>& offsets,
                                                        std::list<uint>& cellstype,
                                                        uint&            firstpoint_globalnumber,
                                                        uint&            last_offset) const
{
    for(const auto& face : kDodecahedronFaces)
    {
        for(int idx : face)
            connectivity.push_back(firstpoint_globalnumber + static_cast<uint>(idx));
        last_offset += 5;
        offsets.push_back(last_offset);
        cellstype.push_back(7);
    }
    firstpoint_globalnumber += 20;
}

// -------------------------------------------------------------------------------------------------
// Returns the shape parameters
template <typename T>
__HOSTDEVICE__ void Dodecahedron<T>::getShapeParameters(T (&params)[5]) const
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
__HOST__ std::string Dodecahedron<T>::getConvexName() const
{
    return "Dodecahedron";
}

// -------------------------------------------------------------------------------------------------
// Writes the dodecahedron in OBJ format
template <typename T>
__HOST__ void Dodecahedron<T>::writeOBJ(std::ostream& f, size_t& firstpoint_number) const
{
    for(int i = 0; i < 20; ++i)
    {
        const Vector3<T> vertex = dodecahedronVertex<T>(i, m_radius);
        f << "v " << vertex[X] << " " << vertex[Y] << " " << vertex[Z] << "\n";
    }

    for(const auto& face : kDodecahedronFaces)
    {
        f << "f";
        for(int idx : face)
            f << " " << firstpoint_number + static_cast<size_t>(idx);
        f << "\n";
    }

    firstpoint_number += 20;
}

// -------------------------------------------------------------------------------------------------
// Explicit template instantiations
template class Dodecahedron<float>;
template class Dodecahedron<double>;