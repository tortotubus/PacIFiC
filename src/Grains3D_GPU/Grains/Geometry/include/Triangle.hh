#ifndef _TRIANGLE_HH_
#define _TRIANGLE_HH_

#include "Convex.hh"
#include "ReaderXML.hh"

// =================================================================================================
/** @brief The class Triangle.

    Convex with the shape of an isosceles triangle lying in the XY plane (Z = 0).
    The base edge (full width = 2*m_halfBase) lies at y = -m_halfHeight and the apex
    is at (0, +m_halfHeight, 0).

    @author A.Yazdani - 2026 - Construction */
// =================================================================================================
template <typename T>
class Triangle : public Convex<T>
{
protected:
    /** @name Parameters */
    //@{
    T m_halfBase;   /**< half-width of the base edge (y = -m_halfHeight) */
    T m_halfHeight; /**< half-height (distance from centre to base/apex) */
    //@}

public:
    /** @name Constructors */
    //@{
    /** @brief Constructor with full base width and full height
        @param base full base width
        @param h    full height */
    __HOSTDEVICE__
    Triangle(T base = T(0), T h = T(0));

    /** @brief Constructor with an input stream
        @param fileIn input stream */
    __HOST__
    Triangle(std::istream& fileIn);

    /** @brief Constructor with an XML node as an input parameter
        @param root XML node (attributes: Base, Height) */
    __HOST__
    Triangle(DOMNode* root);

    /** @brief Destructor */
    __HOSTDEVICE__
    ~Triangle();
    //@}

    /** @name Get methods */
    //@{
    /** @brief Returns the convex type */
    __HOSTDEVICE__
    ConvexType getConvexType() const final;

    /** @brief Returns the half-width of the base edge */
    __HOSTDEVICE__
    T getHalfBase() const;

    /** @brief Returns the half-height */
    __HOSTDEVICE__
    T getHalfHeight() const;
    //@}

    /** @name Methods */
    //@{
    /** @brief Returns a clone of the triangle */
    __HOSTDEVICE__
    Convex<T>* clone() const final;

    /** @brief Returns the triangle area used as geometric measure */
    __HOSTDEVICE__
    T computeVolume() const final;

    /** @brief Computes the inertia tensor (thin planar approximation)
        @param inertia inertia tensor */
    __HOSTDEVICE__
    void computeInertia(T (&inertia)[3]) const final;

    /** @brief Returns the circumscribed radius of the triangle */
    __HOSTDEVICE__
    T computeCircumscribedRadius() const final;

    /** @brief Returns the half-extents of the bounding box fitted to the triangle */
    __HOSTDEVICE__
    Vector3<T> computeBoundingBox() const final;

    /** @brief Returns the tightest bounding cylinder fitted to the triangle in body-local frame */
    __HOSTDEVICE__
    Vector3<T> computeBoundingCylinder() const final;

    /** @brief Triangle support function: returns the support point P satisfying max(P.v)
        @param v direction */
    __HOSTDEVICE__
    Vector3<T> support(const Vector3<T>& v) const final;

    /** @brief Eroded triangle support: evaluates support with uniformly shrunk dimensions
        @param v     direction
        @param crust crust thickness */
    __HOSTDEVICE__
    Vector3<T> support(const Vector3<T>& v, T crust) const final;

    /** @brief Returns whether point p lies in the triangle
        @param p point */
    __HOSTDEVICE__
    bool isInside(const Vector3<T>& p) const final;
    //@}

    /** @name I/O methods */
    //@{
    /** @brief Input operator
        @param fileIn input stream */
    __HOST__
    void readConvex(std::istream& fileIn) final;

    /** @brief Output operator
        @param fileOut output stream */
    __HOST__
    void writeConvex(std::ostream& fileOut) const final;

    /** @brief Returns the number of points to write the triangle in a Paraview format */
    __HOST__
    int numberOfPoints_PARAVIEW() const final;

    /** @brief Returns the number of elementary polytopes to write the triangle in a Paraview
        format */
    __HOST__
    int numberOfCells_PARAVIEW() const final;

    /** @brief Returns a list of points describing the triangle in a Paraview format
        @param transform geometric transformation
        @param translation additional center of mass translation */
    __HOST__
    std::list<Vector3<T>> writePoints_PARAVIEW(const Transform3<T>& transform,
                                               Vector3<T> const*    translation) const final;

    /** @brief Writes the connectivity of the triangle in a Paraview format
        @param connectivity connectivity of Paraview polytopes
        @param offsets connectivity offsets
        @param cellstype Paraview polytopes type
        @param firstpoint_globalnumber global number of the 1st point
        @param last_offset last offset used for the previous convex shape */
    __HOST__
    void writeConnection_PARAVIEW(std::list<uint>& connectivity,
                                  std::list<uint>& offsets,
                                  std::list<uint>& cellstype,
                                  uint&            firstpoint_globalnumber,
                                  uint&            last_offset) const final;

    /** @brief Returns triangle shape parameters: [halfBase, halfHeight, 0, 0, 0] */
    __HOSTDEVICE__
    void getShapeParameters(T (&params)[5]) const final;

    /** @brief Returns "Triangle" */
    __HOST__
    std::string getConvexName() const final;

    /** @brief Triangle has no volumetric OBJ mesh; does nothing */
    __HOST__
    void writeOBJ(std::ostream&, size_t&) const final;
    //@}
};

#endif
