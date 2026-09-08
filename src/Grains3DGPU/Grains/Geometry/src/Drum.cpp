#include "Drum.hh"
#include "VectorMath.hh"

// multiple of 4
#define visuNodeNbOnPer 32

// -------------------------------------------------------------------------------------------------
// Constructor with radius and height as input parameters
template <typename T>
__HOSTDEVICE__ Drum<T>::Drum(T r, T h)
    : m_radius(r)
    , m_halfHeight(h / T(2))
{
}

// -------------------------------------------------------------------------------------------------
// Constructor with an input stream
template <typename T>
__HOST__ Drum<T>::Drum(std::istream& fileIn)
{
    readConvex(fileIn);
}

// -------------------------------------------------------------------------------------------------
// Constructor with an XML node as an input parameter
template <typename T>
__HOST__ Drum<T>::Drum(DOMNode* root)
{
    m_radius     = T(ReaderXML::getNodeAttr_Double(root, "Radius"));
    m_halfHeight = T(ReaderXML::getNodeAttr_Double(root, "Height")) / T(2);
}

// -------------------------------------------------------------------------------------------------
// Destructor
template <typename T>
__HOSTDEVICE__ Drum<T>::~Drum()
{
}

// -------------------------------------------------------------------------------------------------
// Returns the convex type
template <typename T>
__HOSTDEVICE__ ConvexType Drum<T>::getConvexType() const
{
    return (DRUM);
}

// -------------------------------------------------------------------------------------------------
// Returns the radius
template <typename T>
__HOSTDEVICE__ T Drum<T>::getRadius() const
{
    return (m_radius);
}

// -------------------------------------------------------------------------------------------------
// Returns the height
template <typename T>
__HOSTDEVICE__ T Drum<T>::getHeight() const
{
    return (T(2) * m_halfHeight);
}

// -------------------------------------------------------------------------------------------------
// Returns a clone of the drum
template <typename T>
__HOSTDEVICE__ Convex<T>* Drum<T>::clone() const
{
    return (new Drum<T>(m_radius, T(2) * m_halfHeight));
}

// -------------------------------------------------------------------------------------------------
// Returns the volume of the drum
template <typename T>
__HOSTDEVICE__ T Drum<T>::computeVolume() const
{
    return (TWO_PI<T> * m_radius * (T(2) * m_halfHeight));
}

// -------------------------------------------------------------------------------------------------
// Computes the diagonal inertia tensor
template <typename T>
__HOSTDEVICE__ void Drum<T>::computeInertia(T (&inertia)[3]) const
{
    const T height = T(2) * m_halfHeight;
    const T area   = TWO_PI<T> * m_radius * height;
    inertia[0] = inertia[2] = area * (m_radius * m_radius / T(2) + height * height / T(12));
    inertia[1]              = area * m_radius * m_radius;
}

// -------------------------------------------------------------------------------------------------
// Returns the circumscribed radius of the drum
template <typename T>
__HOSTDEVICE__ T Drum<T>::computeCircumscribedRadius() const
{
    return (sqrt(m_radius * m_radius + m_halfHeight * m_halfHeight));
}

// -------------------------------------------------------------------------------------------------
// Returns the bounding box of the drum
template <typename T>
__HOSTDEVICE__ Vector3<T> Drum<T>::computeBoundingBox() const
{
    return (Vector3<T>(m_radius, m_halfHeight, m_radius));
}

// -------------------------------------------------------------------------------------------------
// Returns the bounding cylinder of the drum (exact fit)
template <typename T>
__HOSTDEVICE__ Vector3<T> Drum<T>::computeBoundingCylinder() const
{
    return (Vector3<T>(m_radius, m_halfHeight, T(1)));
}

// -------------------------------------------------------------------------------------------------
// Drum support function, returns the support point P, i.e. the point on
// the surface of the matching finite cylinder that satisfies max(P.v)
template <typename T>
__HOSTDEVICE__ Vector3<T> Drum<T>::support(const Vector3<T>& v) const
{
    const T s  = sqrt(v[X] * v[X] + v[Z] * v[Z]);
    const T hy = fabs(v[Y]) < EPS<T> ? T(0) : (v[Y] < T(0) ? -m_halfHeight : m_halfHeight);
    if(s > EPS<T>)
    {
        const T d = m_radius / s;
        return (Vector3<T>(v[X] * d, hy, v[Z] * d));
    }
    return (Vector3<T>(T(0), hy, T(0)));
}

// -------------------------------------------------------------------------------------------------
// Eroded drum support: evaluates support of matching cylinder with shrunk dimensions
template <typename T>
__HOSTDEVICE__ Vector3<T> Drum<T>::support(const Vector3<T>& v, T crust) const
{
    const T radius     = m_radius - crust;
    const T halfHeight = m_halfHeight - crust;
    const T s          = sqrt(v[X] * v[X] + v[Z] * v[Z]);
    const T hy         = fabs(v[Y]) < EPS<T> ? T(0) : (v[Y] < T(0) ? -halfHeight : halfHeight);
    if(s > EPS<T>)
    {
        const T d = radius / s;
        return (Vector3<T>(v[X] * d, hy, v[Z] * d));
    }
    return (Vector3<T>(T(0), hy, T(0)));
}

// -------------------------------------------------------------------------------------------------
// Input operator
template <typename T>
__HOST__ void Drum<T>::readConvex(std::istream& fileIn)
{
    fileIn >> m_radius >> m_halfHeight;
    m_halfHeight /= T(2);
}

// -------------------------------------------------------------------------------------------------
// Output operator
template <typename T>
__HOST__ void Drum<T>::writeConvex(std::ostream& fileOut) const
{
    fileOut << "Drum: " << m_radius << ", " << T(2) * m_halfHeight << ".\n";
}

// -------------------------------------------------------------------------------------------------
// Returns the number of points to write the drum in a Paraview format
template <typename T>
__HOST__ int Drum<T>::numberOfPoints_PARAVIEW() const
{
    return (2 * visuNodeNbOnPer);
}

// -------------------------------------------------------------------------------------------------
// Returns the number of elementary polytopes to write the drum in a Paraview format
template <typename T>
__HOST__ int Drum<T>::numberOfCells_PARAVIEW() const
{
    return (visuNodeNbOnPer);
}

// -------------------------------------------------------------------------------------------------
// Returns a list of points describing the drum in a Paraview format
template <typename T>
__HOST__ std::list<Vector3<T>> Drum<T>::writePoints_PARAVIEW(const Transform3<T>& transform,
                                                             Vector3<T> const*    translation) const
{
    std::list<Vector3<T>> points;
    Vector3<T>            point;
    const T               dtheta = TWO_PI<T> / T(visuNodeNbOnPer);

    point[Y] = -m_halfHeight;
    for(int i = 0; i < visuNodeNbOnPer; ++i)
    {
        point[X]       = m_radius * cos(T(i) * dtheta);
        point[Z]       = m_radius * sin(T(i) * dtheta);
        Vector3<T> out = transform(point);
        if(translation)
            out += *translation;
        points.push_back(out);
    }

    point[Y] = m_halfHeight;
    for(int i = 0; i < visuNodeNbOnPer; ++i)
    {
        point[X]       = m_radius * cos(T(i) * dtheta);
        point[Z]       = m_radius * sin(T(i) * dtheta);
        Vector3<T> out = transform(point);
        if(translation)
            out += *translation;
        points.push_back(out);
    }

    return (points);
}

// -------------------------------------------------------------------------------------------------
// Writes the connectivity of the drum in a Paraview format
template <typename T>
__HOST__ void Drum<T>::writeConnection_PARAVIEW(std::list<uint>& connectivity,
                                                std::list<uint>& offsets,
                                                std::list<uint>& cellstype,
                                                uint&            firstpoint_globalnumber,
                                                uint&            last_offset) const
{
    for(int i = 0; i < visuNodeNbOnPer; ++i)
    {
        const uint lower0 = firstpoint_globalnumber + static_cast<uint>(i);
        const uint lower1 = firstpoint_globalnumber + static_cast<uint>((i + 1) % visuNodeNbOnPer);
        const uint upper1 = firstpoint_globalnumber + static_cast<uint>(visuNodeNbOnPer)
                            + static_cast<uint>((i + 1) % visuNodeNbOnPer);
        const uint upper0
            = firstpoint_globalnumber + static_cast<uint>(visuNodeNbOnPer) + static_cast<uint>(i);
        connectivity.push_back(lower0);
        connectivity.push_back(lower1);
        connectivity.push_back(upper1);
        connectivity.push_back(upper0);
        last_offset += 4;
        offsets.push_back(last_offset);
        cellstype.push_back(9);
    }

    firstpoint_globalnumber += static_cast<uint>(2 * visuNodeNbOnPer);
}

// -------------------------------------------------------------------------------------------------
// Returns drum shape parameters: [radius, halfHeight, 0, 0, 0]
template <typename T>
__HOSTDEVICE__ void Drum<T>::getShapeParameters(T (&params)[5]) const
{
    params[0] = m_radius;
    params[1] = m_halfHeight;
    params[2] = T(0);
    params[3] = T(0);
    params[4] = T(0);
}

// -------------------------------------------------------------------------------------------------
// Returns the convex name
template <typename T>
__HOST__ std::string Drum<T>::getConvexName() const
{
    return ("Drum");
}

// -------------------------------------------------------------------------------------------------
// Writes the drum as a Wavefront OBJ mesh
template <typename T>
__HOST__ void Drum<T>::writeOBJ(std::ostream& f, size_t& firstpoint_number) const
{
    const T dtheta = TWO_PI<T> / T(visuNodeNbOnPer);
    for(int ring = 0; ring < 2; ++ring)
    {
        const T y = ring == 0 ? -m_halfHeight : m_halfHeight;
        for(int i = 0; i < visuNodeNbOnPer; ++i)
        {
            f << "v " << m_radius * cos(T(i) * dtheta) << " " << y << " "
              << m_radius * sin(T(i) * dtheta) << "\n";
        }
    }

    for(int i = 0; i < visuNodeNbOnPer; ++i)
    {
        const size_t lower0 = firstpoint_number + static_cast<size_t>(i);
        const size_t lower1 = firstpoint_number + static_cast<size_t>((i + 1) % visuNodeNbOnPer);
        const size_t upper1 = firstpoint_number + static_cast<size_t>(visuNodeNbOnPer)
                              + static_cast<size_t>((i + 1) % visuNodeNbOnPer);
        const size_t upper0
            = firstpoint_number + static_cast<size_t>(visuNodeNbOnPer) + static_cast<size_t>(i);
        f << "f " << lower0 << " " << lower1 << " " << upper1 << " " << upper0 << "\n";
    }

    firstpoint_number += static_cast<size_t>(2 * visuNodeNbOnPer);
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class Drum<float>;
template class Drum<double>;