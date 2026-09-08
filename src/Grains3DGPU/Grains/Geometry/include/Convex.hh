#ifndef _CONVEX_HH_
#define _CONVEX_HH_

#include "Transform3.hh"
#include <string>

// Convex types
enum ConvexType
{
    SPHERE       = 0,
    BOX          = 1,
    OCTAHEDRON   = 2,
    DODECAHEDRON = 3,
    CYLINDER     = 4,
    CONE         = 5,
    SUPERQUADRIC = 6,
    RECTANGLE    = 7,
    POLYHEDRON   = 8,
    TRAPEZOID    = 9,
    TRIANGLE     = 10,
    DRUM         = 11
};

// =================================================================================================
/** @brief The class Convex. Convex bodies - The base class for various particle shapes.

    @author A.Yazdani - 2023 - Construction
    @author A.Yazdani - 2024 - Modificiation */
// =================================================================================================
template <typename T>
class Convex
{
protected:
    /**@name Contructors */
    //@{
    /** @brief Default constructor (forbidden except in derived classes) */
    __HOSTDEVICE__
    Convex();
    //@}

public:
    /** @name Constructors */
    //@{
    /** @brief Destructor */
    __HOSTDEVICE__
    virtual ~Convex();
    //@}

    /** @name Get methods */
    //@{
    /** @brief Returns the convex type */
    __HOSTDEVICE__
    virtual ConvexType getConvexType() const = 0;

    /** @brief Returns a human-readable name for the convex type (e.g. "Sphere", "Box") */
    __HOST__
    virtual std::string getConvexName() const = 0;
    //@}

    /** @name Methods */
    //@{
    /** @brief Returns a clone of the convex */
    __HOSTDEVICE__
    virtual Convex<T>* clone() const = 0;

    /** @brief Returns the volume of the convex shape */
    __HOSTDEVICE__
    virtual T computeVolume() const = 0;

    /** @brief Computes the diagonal inertia tensor
        @param inertia diagonal inertia tensor (3 components: Ixx, Iyy, Izz) */
    __HOSTDEVICE__
    virtual void computeInertia(T (&inertia)[3]) const = 0;

    /** @brief Computes and returns the circumscribed radius of the reference convex shape */
    __HOSTDEVICE__
    virtual T computeCircumscribedRadius() const = 0;

    /** @brief Returns the half-length of the bounding box fitted to the convex without considering
        the transformation */
    __HOSTDEVICE__
    virtual Vector3<T> computeBoundingBox() const = 0;

    /** @brief Returns the tightest bounding cylinder fitted to the convex in body-local frame.
        The returned Vector3 encodes: [0] radius, [1] half-height, [2] axis index (0=X, 1=Y, 2=Z).
     */
    __HOSTDEVICE__
    virtual Vector3<T> computeBoundingCylinder() const = 0;

    /** @brief Convex support function, returns the support point P, i.e. the point on the surface
        of the convex shape that satisfies max(P.v)
        @param v direction vector */
    __HOSTDEVICE__
    virtual Vector3<T> support(const Vector3<T>& v) const = 0;

    /** @brief Eroded convex support function, returns the support point of the
        shape shrunk inward by crust. The default implementation applies generic Minkowski erosion:
        support(v) - crust/norm(v) * v. Built-in shape subclasses override this to evaluate the
        support of the dimensionally-reduced shape directly, preserving the exact shape geometry.
        @param v direction vector
        @param crust erosion thickness */
    __HOSTDEVICE__
    virtual Vector3<T> support(const Vector3<T>& v, T crust) const;

    /** @brief Returns whether point p lies in the convex shape
        @param p point */
    __HOSTDEVICE__
    virtual bool isInside(const Vector3<T>& p) const;
    //@}

    /** @name I/O methods */
    //@{
    /** @brief Input operator
        @param fileIn input stream */
    __HOST__
    virtual void readConvex(std::istream& fileIn) = 0;

    /** @brief Output operator
        @param fileOut output stream */
    __HOST__
    virtual void writeConvex(std::ostream& fileOut) const = 0;

    /** @brief Returns the number of points to write the convex in a Paraview format */
    __HOST__
    virtual int numberOfPoints_PARAVIEW() const = 0;

    /** @brief Returns the number of elementary polytopes to write the convex in a Paraview format
     */
    __HOST__
    virtual int numberOfCells_PARAVIEW() const = 0;

    /** @brief Writes the list of points describing the convex to an stream
        @param f output stream
        @param transform geometric transformation
        @param translation additional center of mass translation */
    __HOST__
    void writePoints_PARAVIEW(std::ostream&        f,
                              const Transform3<T>& transform,
                              const Vector3<T>*    translation = NULL) const;

    /** @brief Returns a list of points describing the convex in a Paraview format
        @param transform geometric transformation
        @param translation additional center of mass translation */
    __HOST__
    virtual std::list<Vector3<T>> writePoints_PARAVIEW(const Transform3<T>& transform,
                                                       Vector3<T> const*    translation) const
        = 0;

    /** @brief Writes the connectivity of the convex in a Paraview format
        @param connectivity connectivity of Paraview polytopes
        @param offsets connectivity offsets
        @param cellstype Paraview polytopes type
        @param firstpoint_globalnumber global number of the 1st point
        @param last_offset last offset used for the previous convex shape */
    __HOST__
    virtual void writeConnection_PARAVIEW(std::list<uint>& connectivity,
                                          std::list<uint>& offsets,
                                          std::list<uint>& cellstype,
                                          uint&            firstpoint_globalnumber,
                                          uint&            last_offset) const
        = 0;

    /** @brief Returns shape-specific geometric parameters.
            Fills a fixed-size array with shape-specific parameters. Default implementation fills
    zeros. Used by lightweight post-processing and ShapeData packing. Layout depends on
    getConvexType():
        Sphere:       [radius, 0, 0, 0, 0]
        Box:          [Lx, Ly, Lz, 0, 0]
        Octahedron:   [radius, 0, 0, 0, 0]
        Dodecahedron: [radius, 0, 0, 0, 0]
        Cylinder:     [radius, height, 0, 0, 0]
        Cone:         [bottomRadius, height, 0, 0, 0]
        Superquadric: [a, b, c, n1, n2]
        Rectangle:    [Lx, Ly, 0, 0, 0]
        Trapezoid:    [halfWidthBottom, halfWidthTop, halfHeight, 0, 0]
        Triangle:     [halfBase, halfHeight, 0, 0, 0]
        Drum:         [radius, halfHeight, 0, 0, 0]
    @param params output array of 5 values */
    __HOSTDEVICE__
    virtual void getShapeParameters(T (&params)[5]) const;

    /** @brief Writes the convex as a Wavefront OBJ mesh to a stream (surface-only, no GC).
        OBJ vertex indices are 1-based. firstpoint_number must be 1 for a single-shape file and
        is advanced by the number of vertices written, matching the reference write_convex_OBJ
        convention so multiple shapes can share a single stream if needed.
        @param f output stream
        @param firstpoint_number 1-based index of the first vertex written; advanced on return */
    __HOST__
    virtual void writeOBJ(std::ostream& f, size_t& firstpoint_number) const = 0;
    //@}
};

/** @name External Methods - I/O methods */
//@{
/** @brief Output operator for Convex: delegates to virtual writeConvex
    @param fileOut output stream
    @param convex convex object */
template <typename T>
__HOST__ std::ostream& operator<<(std::ostream& fileOut, const Convex<T>& convex);

/** @brief Input operator for Convex: delegates to virtual readConvex
    @param fileIn input stream
    @param convex convex object */
template <typename T>
__HOST__ std::istream& operator>>(std::istream& fileIn, Convex<T>& convex);
//@}

#endif