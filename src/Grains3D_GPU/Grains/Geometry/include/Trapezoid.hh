#ifndef _TRAPEZOID_HH_
#define _TRAPEZOID_HH_

#include "Convex.hh"
#include "ReaderXML.hh"

// =================================================================================================
/** @brief The class Trapezoid.

    Convex with the shape of an isosceles trapezoid lying in the XY plane (Z = 0).
    The bottom edge (full width = 2*m_halfWidthBottom) is at y = -m_halfHeight and the top edge
    (full width = 2*m_halfWidthTop) is at y = +m_halfHeight.

    @author A.Yazdani - 2026 - Construction */
// =================================================================================================
template <typename T>
class Trapezoid : public Convex<T>
{
protected:
    /** @name Parameters */
    //@{
    T m_halfWidthBottom; /**< half-width at the bottom edge (y = -m_halfHeight) */
    T m_halfWidthTop;    /**< half-width at the top edge    (y = +m_halfHeight) */
    T m_halfHeight;      /**< half-height of the trapezoid                       */
    //@}

public:
    /** @name Constructors */
    //@{
    /** @brief Constructor with full widths and full height
        @param wb full bottom width
        @param wt full top width
        @param h  full height */
    __HOSTDEVICE__
    Trapezoid(T wb = T(0), T wt = T(0), T h = T(0));

    /** @brief Constructor with an input stream
        @param fileIn input stream */
    __HOST__
    Trapezoid(std::istream& fileIn);

    /** @brief Constructor with an XML node as an input parameter
        @param root XML node (attributes: WidthBottom, WidthTop, Height) */
    __HOST__
    Trapezoid(DOMNode* root);

    /** @brief Destructor */
    __HOSTDEVICE__
    ~Trapezoid();
    //@}

    /** @name Get methods */
    //@{
    /** @brief Returns the convex type */
    __HOSTDEVICE__
    ConvexType getConvexType() const final;

    /** @brief Returns the half-width at the bottom edge */
    __HOSTDEVICE__
    T getHalfWidthBottom() const;

    /** @brief Returns the half-width at the top edge */
    __HOSTDEVICE__
    T getHalfWidthTop() const;

    /** @brief Returns the half-height */
    __HOSTDEVICE__
    T getHalfHeight() const;
    //@}

    /** @name Methods */
    //@{
    /** @brief Returns a clone of the trapezoid */
    __HOSTDEVICE__
    Convex<T>* clone() const final;

    /** @brief Returns the trapezoid area used as geometric measure */
    __HOSTDEVICE__
    T computeVolume() const final;

    /** @brief Computes the inertia tensor (thin planar approximation)
        @param inertia inertia tensor */
    __HOSTDEVICE__
    void computeInertia(T (&inertia)[3]) const final;

    /** @brief Returns the circumscribed radius of the trapezoid */
    __HOSTDEVICE__
    T computeCircumscribedRadius() const final;

    /** @brief Returns the half-extents of the bounding box fitted to the trapezoid */
    __HOSTDEVICE__
    Vector3<T> computeBoundingBox() const final;

    /** @brief Returns the tightest bounding cylinder fitted to the trapezoid in body-local frame */
    __HOSTDEVICE__
    Vector3<T> computeBoundingCylinder() const final;

    /** @brief Trapezoid support function: returns the support point P satisfying max(P.v)
        @param v direction */
    __HOSTDEVICE__
    Vector3<T> support(const Vector3<T>& v) const final;

    /** @brief Eroded trapezoid support: evaluates support with uniformly shrunk dimensions
        @param v     direction
        @param crust crust thickness */
    __HOSTDEVICE__
    Vector3<T> support(const Vector3<T>& v, T crust) const final;

    /** @brief Returns whether point p lies in the trapezoid
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

    /** @brief Returns the number of points to write the trapezoid in a Paraview format */
    __HOST__
    int numberOfPoints_PARAVIEW() const final;

    /** @brief Returns the number of elementary polytopes to write the trapezoid in a Paraview
        format */
    __HOST__
    int numberOfCells_PARAVIEW() const final;

    /** @brief Returns a list of points describing the trapezoid in a Paraview format
        @param transform geometric transformation
        @param translation additional center of mass translation */
    __HOST__
    std::list<Vector3<T>> writePoints_PARAVIEW(const Transform3<T>& transform,
                                               Vector3<T> const*    translation) const final;

    /** @brief Writes the connectivity of the trapezoid in a Paraview format
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

    /** @brief Returns trapezoid shape parameters: [halfWidthBottom, halfWidthTop, halfHeight, 0, 0]
     */
    __HOSTDEVICE__
    void getShapeParameters(T (&params)[5]) const final;

    /** @brief Returns "Trapezoid" */
    __HOST__
    std::string getConvexName() const final;

    /** @brief Trapezoid has no volumetric OBJ mesh; does nothing */
    __HOST__
    void writeOBJ(std::ostream&, size_t&) const final;
    //@}
};

#endif
