#include "Triangle.hh"

// -------------------------------------------------------------------------------------------------
// Constructor with full base width and full height
template <typename T>
__HOSTDEVICE__ Triangle<T>::Triangle(T base, T h)
    : m_halfBase(base / T(2))
    , m_halfHeight(h / T(2))
{
}

// -------------------------------------------------------------------------------------------------
// Constructor with an input stream
template <typename T>
__HOST__ Triangle<T>::Triangle(std::istream& fileIn)
{
    readConvex(fileIn);
}

// -------------------------------------------------------------------------------------------------
// Constructor with an XML node as an input parameter
template <typename T>
__HOST__ Triangle<T>::Triangle(DOMNode* root)
{
    m_halfBase   = T(ReaderXML::getNodeAttr_Double(root, "Base")) / T(2);
    m_halfHeight = T(ReaderXML::getNodeAttr_Double(root, "Height")) / T(2);
}

// -------------------------------------------------------------------------------------------------
// Destructor
template <typename T>
__HOSTDEVICE__ Triangle<T>::~Triangle()
{
}

// -------------------------------------------------------------------------------------------------
// Returns the convex type
template <typename T>
__HOSTDEVICE__ ConvexType Triangle<T>::getConvexType() const
{
    return (TRIANGLE);
}

// -------------------------------------------------------------------------------------------------
// Returns the half-width of the base edge
template <typename T>
__HOSTDEVICE__ T Triangle<T>::getHalfBase() const
{
    return (m_halfBase);
}

// -------------------------------------------------------------------------------------------------
// Returns the half-height
template <typename T>
__HOSTDEVICE__ T Triangle<T>::getHalfHeight() const
{
    return (m_halfHeight);
}

// -------------------------------------------------------------------------------------------------
// Returns a clone of the triangle
template <typename T>
__HOSTDEVICE__ Convex<T>* Triangle<T>::clone() const
{
    return (new Triangle<T>(T(2) * m_halfBase, T(2) * m_halfHeight));
}

// -------------------------------------------------------------------------------------------------
// Returns the area of the triangle: 0.5 * base * height = 0.5 * (2*hb) * (2*hh) = 2*hb*hh
template <typename T>
__HOSTDEVICE__ T Triangle<T>::computeVolume() const
{
    return (T(2) * m_halfBase * m_halfHeight);
}

// -------------------------------------------------------------------------------------------------
// Computes the diagonal inertia tensor (thin planar shape, bounding-rectangle approximation)
template <typename T>
__HOSTDEVICE__ void Triangle<T>::computeInertia(T (&inertia)[3]) const
{
    const T hb = m_halfBase;
    const T hh = m_halfHeight;
    inertia[0] = T(4) / T(3) * hb * hh * hh * hh;
    inertia[1] = T(4) / T(3) * hb * hb * hb * hh;
    inertia[2] = T(4) / T(3) * hb * hh * (hb * hb + hh * hh);
}

// -------------------------------------------------------------------------------------------------
// Returns the circumscribed radius: max distance from origin to any vertex.
// Vertices: (±hb, -hh, 0) and (0, +hh, 0).
template <typename T>
__HOSTDEVICE__ T Triangle<T>::computeCircumscribedRadius() const
{
    const T r2base = m_halfBase * m_halfBase + m_halfHeight * m_halfHeight;
    const T r2apex = m_halfHeight * m_halfHeight;
    return sqrt(r2base > r2apex ? r2base : r2apex);
}

// -------------------------------------------------------------------------------------------------
// Returns the half-extents of the bounding box fitted to the triangle
template <typename T>
__HOSTDEVICE__ Vector3<T> Triangle<T>::computeBoundingBox() const
{
    return (Vector3<T>(m_halfBase, m_halfHeight, EPS<T>));
}

// -------------------------------------------------------------------------------------------------
// Returns the bounding cylinder fitted to the triangle in body-local frame
template <typename T>
__HOSTDEVICE__ Vector3<T> Triangle<T>::computeBoundingCylinder() const
{
    return Vector3<T>(computeCircumscribedRadius(), EPS<T>, T(2));
}

// -------------------------------------------------------------------------------------------------
// Triangle support function. Returns the support point P satisfying max(P.v).
// Vertices: (-hb, -hh, 0), (+hb, -hh, 0), (0, +hh, 0)
// Compare apex score = v[Y]*hh vs best base-corner score = |v[X]|*hb - v[Y]*hh.
// Apex wins when v[Y]*hh > |v[X]|*hb - v[Y]*hh, i.e. 2*v[Y]*hh > |v[X]|*hb.
template <typename T>
__HOSTDEVICE__ Vector3<T> Triangle<T>::support(const Vector3<T>& v) const
{
    const T scoreApex = v[Y] * m_halfHeight;
    const T scoreBase = -v[Y] * m_halfHeight + fabs(v[X]) * m_halfBase;
    if(scoreApex >= scoreBase)
        return (Vector3<T>(T(0), m_halfHeight, T(0)));
    return (Vector3<T>(v[X] < T(0) ? -m_halfBase : m_halfBase, -m_halfHeight, T(0)));
}

// -------------------------------------------------------------------------------------------------
// Eroded triangle support: evaluates support with uniformly shrunk dimensions
template <typename T>
__HOSTDEVICE__ Vector3<T> Triangle<T>::support(const Vector3<T>& v, T crust) const
{
    const T hb_e      = m_halfBase - crust > T(0) ? m_halfBase - crust : T(0);
    const T hh_e      = m_halfHeight - crust > T(0) ? m_halfHeight - crust : T(0);
    const T scoreApex = v[Y] * hh_e;
    const T scoreBase = -v[Y] * hh_e + fabs(v[X]) * hb_e;
    if(scoreApex >= scoreBase)
        return (Vector3<T>(T(0), hh_e, T(0)));
    return (Vector3<T>(v[X] < T(0) ? -hb_e : hb_e, -hh_e, T(0)));
}

// -------------------------------------------------------------------------------------------------
// Returns whether point p lies in the triangle
template <typename T>
__HOSTDEVICE__ bool Triangle<T>::isInside(const Vector3<T>& p) const
{
    if(p[Y] < -m_halfHeight || p[Y] > m_halfHeight)
        return false;
    // Half-width at y: hb * (hh - y) / (2*hh)  (0 at apex, hb at base)
    const T hw = m_halfBase * (m_halfHeight - p[Y]) / (T(2) * m_halfHeight);
    return (p[X] >= -hw && p[X] <= hw);
}

// -------------------------------------------------------------------------------------------------
// Input operator
template <typename T>
__HOST__ void Triangle<T>::readConvex(std::istream& fileIn)
{
    T base, h;
    fileIn >> base >> h;
    m_halfBase   = base / T(2);
    m_halfHeight = h / T(2);
}

// -------------------------------------------------------------------------------------------------
// Output operator
template <typename T>
__HOST__ void Triangle<T>::writeConvex(std::ostream& fileOut) const
{
    fileOut << "Triangle: " << T(2) * m_halfBase << ", " << T(2) * m_halfHeight << ".\n";
}

// -------------------------------------------------------------------------------------------------
// Returns the number of points to write the triangle in a Paraview format
template <typename T>
__HOST__ int Triangle<T>::numberOfPoints_PARAVIEW() const
{
    return (3);
}

// -------------------------------------------------------------------------------------------------
// Returns the number of elementary polytopes to write the triangle in a Paraview format
template <typename T>
__HOST__ int Triangle<T>::numberOfCells_PARAVIEW() const
{
    return (1);
}

// -------------------------------------------------------------------------------------------------
// Returns a list of points describing the triangle in a Paraview format.
// Vertices in order: base-left, base-right, apex
template <typename T>
__HOST__ std::list<Vector3<T>>
         Triangle<T>::writePoints_PARAVIEW(const Transform3<T>& transform,
                                      Vector3<T> const*    translation) const
{
    std::list<Vector3<T>> ParaviewPoints;
    Vector3<T>            p;
    p.setValue(-m_halfBase, -m_halfHeight, T(0));
    ParaviewPoints.push_back(transform(p));
    p.setValue(m_halfBase, -m_halfHeight, T(0));
    ParaviewPoints.push_back(transform(p));
    p.setValue(T(0), m_halfHeight, T(0));
    ParaviewPoints.push_back(transform(p));
    return (ParaviewPoints);
}

// -------------------------------------------------------------------------------------------------
// Writes the connectivity of the triangle in a Paraview format (triangle cell type 5)
template <typename T>
__HOST__ void Triangle<T>::writeConnection_PARAVIEW(std::list<uint>& connectivity,
                                                    std::list<uint>& offsets,
                                                    std::list<uint>& cellstype,
                                                    uint&            firstpoint_globalnumber,
                                                    uint&            last_offset) const
{
    uint count = firstpoint_globalnumber;
    for(uint i = 0; i < 3; ++i)
    {
        connectivity.push_back(count);
        ++count;
    }
    last_offset += 3;
    offsets.push_back(last_offset);
    cellstype.push_back(5);  // VTK_TRIANGLE

    firstpoint_globalnumber += 3;
}

// -------------------------------------------------------------------------------------------------
// Returns triangle shape parameters: [halfBase, halfHeight, 0, 0, 0]
template <typename T>
__HOSTDEVICE__ void Triangle<T>::getShapeParameters(T (&params)[5]) const
{
    params[0] = m_halfBase;
    params[1] = m_halfHeight;
    params[2] = T(0);
    params[3] = T(0);
    params[4] = T(0);
}

// -------------------------------------------------------------------------------------------------
// Returns the convex name
template <typename T>
__HOST__ std::string Triangle<T>::getConvexName() const
{
    return "Triangle";
}

// -------------------------------------------------------------------------------------------------
// OBJ mesh: Triangle is a 2D shape; nothing to write
template <typename T>
__HOST__ void Triangle<T>::writeOBJ(std::ostream&, size_t&) const
{
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class Triangle<float>;
template class Triangle<double>;
