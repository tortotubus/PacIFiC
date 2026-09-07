#include "Trapezoid.hh"

// -------------------------------------------------------------------------------------------------
// Constructor with full widths and full height
template <typename T>
__HOSTDEVICE__ Trapezoid<T>::Trapezoid(T wb, T wt, T h)
    : m_halfWidthBottom(wb / T(2))
    , m_halfWidthTop(wt / T(2))
    , m_halfHeight(h / T(2))
{
}

// -------------------------------------------------------------------------------------------------
// Constructor with an input stream
template <typename T>
__HOST__ Trapezoid<T>::Trapezoid(std::istream& fileIn)
{
    readConvex(fileIn);
}

// -------------------------------------------------------------------------------------------------
// Constructor with an XML node as an input parameter
template <typename T>
__HOST__ Trapezoid<T>::Trapezoid(DOMNode* root)
{
    m_halfWidthBottom = T(ReaderXML::getNodeAttr_Double(root, "WidthBottom")) / T(2);
    m_halfWidthTop    = T(ReaderXML::getNodeAttr_Double(root, "WidthTop")) / T(2);
    m_halfHeight      = T(ReaderXML::getNodeAttr_Double(root, "Height")) / T(2);
}

// -------------------------------------------------------------------------------------------------
// Destructor
template <typename T>
__HOSTDEVICE__ Trapezoid<T>::~Trapezoid()
{
}

// -------------------------------------------------------------------------------------------------
// Returns the convex type
template <typename T>
__HOSTDEVICE__ ConvexType Trapezoid<T>::getConvexType() const
{
    return (TRAPEZOID);
}

// -------------------------------------------------------------------------------------------------
// Returns the half-width at the bottom edge
template <typename T>
__HOSTDEVICE__ T Trapezoid<T>::getHalfWidthBottom() const
{
    return (m_halfWidthBottom);
}

// -------------------------------------------------------------------------------------------------
// Returns the half-width at the top edge
template <typename T>
__HOSTDEVICE__ T Trapezoid<T>::getHalfWidthTop() const
{
    return (m_halfWidthTop);
}

// -------------------------------------------------------------------------------------------------
// Returns the half-height
template <typename T>
__HOSTDEVICE__ T Trapezoid<T>::getHalfHeight() const
{
    return (m_halfHeight);
}

// -------------------------------------------------------------------------------------------------
// Returns a clone of the trapezoid
template <typename T>
__HOSTDEVICE__ Convex<T>* Trapezoid<T>::clone() const
{
    return (new Trapezoid<T>(T(2) * m_halfWidthBottom, T(2) * m_halfWidthTop, T(2) * m_halfHeight));
}

// -------------------------------------------------------------------------------------------------
// Returns the area of the trapezoid: (halfWidthBottom + halfWidthTop) * 2 * (2 * halfHeight) / 2
// = (wb + wt) * 2 * H where wb, wt are half-widths and H = halfHeight
// i.e., area = (2wb + 2wt) / 2 * 2H = (wb + wt) * 2H
template <typename T>
__HOSTDEVICE__ T Trapezoid<T>::computeVolume() const
{
    return ((m_halfWidthBottom + m_halfWidthTop) * T(2) * m_halfHeight);
}

// -------------------------------------------------------------------------------------------------
// Computes the diagonal inertia tensor (thin planar shape in XY, rotation about Z dominates)
template <typename T>
__HOSTDEVICE__ void Trapezoid<T>::computeInertia(T (&inertia)[3]) const
{
    // Use bounding-rectangle approximation for inertia
    const T wb   = m_halfWidthBottom;
    const T wt   = m_halfWidthTop;
    const T H    = m_halfHeight;
    const T Wmax = wb > wt ? wb : wt;
    // Approximate as bounding rectangle [2Wmax x 2H]
    inertia[0] = T(4) / T(3) * Wmax * H * H * H;
    inertia[1] = T(4) / T(3) * Wmax * Wmax * Wmax * H;
    inertia[2] = T(4) / T(3) * Wmax * H * (Wmax * Wmax + H * H);
}

// -------------------------------------------------------------------------------------------------
// Returns the circumscribed radius of the trapezoid
template <typename T>
__HOSTDEVICE__ T Trapezoid<T>::computeCircumscribedRadius() const
{
    // Furthest vertex from origin: max(|bottom-corner|, |top-corner|)
    const T r2bot = m_halfWidthBottom * m_halfWidthBottom + m_halfHeight * m_halfHeight;
    const T r2top = m_halfWidthTop * m_halfWidthTop + m_halfHeight * m_halfHeight;
    return sqrt(r2bot > r2top ? r2bot : r2top);
}

// -------------------------------------------------------------------------------------------------
// Returns the half-extents of the bounding box fitted to the trapezoid
template <typename T>
__HOSTDEVICE__ Vector3<T> Trapezoid<T>::computeBoundingBox() const
{
    const T maxW = m_halfWidthBottom > m_halfWidthTop ? m_halfWidthBottom : m_halfWidthTop;
    return (Vector3<T>(maxW, m_halfHeight, EPS<T>));
}

// -------------------------------------------------------------------------------------------------
// Returns the bounding cylinder fitted to the trapezoid in body-local frame
template <typename T>
__HOSTDEVICE__ Vector3<T> Trapezoid<T>::computeBoundingCylinder() const
{
    const T radius = computeCircumscribedRadius();
    return Vector3<T>(radius, EPS<T>, T(2));
}

// -------------------------------------------------------------------------------------------------
// Trapezoid support function. Returns the support point P satisfying max(P.v).
// The four vertices are: (-wb, -hh, 0), (+wb, -hh, 0), (-wt, +hh, 0), (+wt, +hh, 0)
// We pick the face (top or bottom) that gives the higher score, then the ±X vertex on that face.
template <typename T>
__HOSTDEVICE__ Vector3<T> Trapezoid<T>::support(const Vector3<T>& v) const
{
    // Score for each face: dot(v, vertex_on_face)
    // Top face vertex in v-direction: v[X]*(±wt) + v[Y]*hh -> best X sign gives |v[X]|*wt
    // Bottom face vertex: |v[X]|*wb - v[Y]*hh
    const T scoreTop = v[Y] * m_halfHeight + fabs(v[X]) * m_halfWidthTop;
    const T scorBot  = -v[Y] * m_halfHeight + fabs(v[X]) * m_halfWidthBottom;
    const T hy       = (scoreTop >= scorBot) ? m_halfHeight : -m_halfHeight;
    const T hw       = (scoreTop >= scorBot) ? m_halfWidthTop : m_halfWidthBottom;
    return (Vector3<T>(v[X] < T(0) ? -hw : hw, hy, T(0)));
}

// -------------------------------------------------------------------------------------------------
// Eroded trapezoid support: evaluates support with uniformly shrunk dimensions
template <typename T>
__HOSTDEVICE__ Vector3<T> Trapezoid<T>::support(const Vector3<T>& v, T crust) const
{
    const T wb_e     = m_halfWidthBottom - crust > T(0) ? m_halfWidthBottom - crust : T(0);
    const T wt_e     = m_halfWidthTop - crust > T(0) ? m_halfWidthTop - crust : T(0);
    const T hh_e     = m_halfHeight - crust > T(0) ? m_halfHeight - crust : T(0);
    const T scoreTop = v[Y] * hh_e + fabs(v[X]) * wt_e;
    const T scorBot  = -v[Y] * hh_e + fabs(v[X]) * wb_e;
    const T hy       = (scoreTop >= scorBot) ? hh_e : -hh_e;
    const T hw       = (scoreTop >= scorBot) ? wt_e : wb_e;
    return (Vector3<T>(v[X] < T(0) ? -hw : hw, hy, T(0)));
}

// -------------------------------------------------------------------------------------------------
// Returns whether point p lies in the trapezoid
template <typename T>
__HOSTDEVICE__ bool Trapezoid<T>::isInside(const Vector3<T>& p) const
{
    if(p[Y] < -m_halfHeight || p[Y] > m_halfHeight)
        return false;
    // Interpolated half-width at y: wb + (wt - wb) * (y + hh) / (2*hh)
    const T t = (p[Y] + m_halfHeight) / (T(2) * m_halfHeight);
    const T w = m_halfWidthBottom + (m_halfWidthTop - m_halfWidthBottom) * t;
    return (p[X] >= -w && p[X] <= w);
}

// -------------------------------------------------------------------------------------------------
// Input operator
template <typename T>
__HOST__ void Trapezoid<T>::readConvex(std::istream& fileIn)
{
    T wb, wt, h;
    fileIn >> wb >> wt >> h;
    m_halfWidthBottom = wb / T(2);
    m_halfWidthTop    = wt / T(2);
    m_halfHeight      = h / T(2);
}

// -------------------------------------------------------------------------------------------------
// Output operator
template <typename T>
__HOST__ void Trapezoid<T>::writeConvex(std::ostream& fileOut) const
{
    fileOut << "Trapezoid: " << T(2) * m_halfWidthBottom << ", " << T(2) * m_halfWidthTop << ", "
            << T(2) * m_halfHeight << ".\n";
}

// -------------------------------------------------------------------------------------------------
// Returns the number of points to write the trapezoid in a Paraview format
template <typename T>
__HOST__ int Trapezoid<T>::numberOfPoints_PARAVIEW() const
{
    return (4);
}

// -------------------------------------------------------------------------------------------------
// Returns the number of elementary polytopes to write the trapezoid in a Paraview format
template <typename T>
__HOST__ int Trapezoid<T>::numberOfCells_PARAVIEW() const
{
    return (1);
}

// -------------------------------------------------------------------------------------------------
// Returns a list of points describing the trapezoid in a Paraview format.
// Vertices in order: bottom-left, bottom-right, top-right, top-left
template <typename T>
__HOST__ std::list<Vector3<T>>
         Trapezoid<T>::writePoints_PARAVIEW(const Transform3<T>& transform,
                                       Vector3<T> const*    translation) const
{
    std::list<Vector3<T>> ParaviewPoints;
    Vector3<T>            p;
    p.setValue(-m_halfWidthBottom, -m_halfHeight, T(0));
    ParaviewPoints.push_back(transform(p));
    p.setValue(m_halfWidthBottom, -m_halfHeight, T(0));
    ParaviewPoints.push_back(transform(p));
    p.setValue(m_halfWidthTop, m_halfHeight, T(0));
    ParaviewPoints.push_back(transform(p));
    p.setValue(-m_halfWidthTop, m_halfHeight, T(0));
    ParaviewPoints.push_back(transform(p));
    return (ParaviewPoints);
}

// -------------------------------------------------------------------------------------------------
// Writes the connectivity of the trapezoid in a Paraview format (quad cell type 8)
template <typename T>
__HOST__ void Trapezoid<T>::writeConnection_PARAVIEW(std::list<uint>& connectivity,
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
// Returns trapezoid shape parameters: [halfWidthBottom, halfWidthTop, halfHeight, 0, 0]
template <typename T>
__HOSTDEVICE__ void Trapezoid<T>::getShapeParameters(T (&params)[5]) const
{
    params[0] = m_halfWidthBottom;
    params[1] = m_halfWidthTop;
    params[2] = m_halfHeight;
    params[3] = T(0);
    params[4] = T(0);
}

// -------------------------------------------------------------------------------------------------
// Returns the convex name
template <typename T>
__HOST__ std::string Trapezoid<T>::getConvexName() const
{
    return "Trapezoid";
}

// -------------------------------------------------------------------------------------------------
// OBJ mesh: Trapezoid is a 2D shape; nothing to write
template <typename T>
__HOST__ void Trapezoid<T>::writeOBJ(std::ostream&, size_t&) const
{
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class Trapezoid<float>;
template class Trapezoid<double>;
