#ifndef _DRUM_HH_
#define _DRUM_HH_

#include "Convex.hh"
#include "ReaderXML.hh"

// =================================================================================================
/** @brief The class Drum.

    Convex proxy for a hollow cylindrical drum.
    Physical contact on the inside of the drum is handled analytically in CollisionDetection.

    @author A.Yazdani - 2026 - Construction */
// =================================================================================================
template <typename T>
class Drum : public Convex<T>
{
protected:
    /** @name Parameters */
    //@{
    T m_radius;
    T m_halfHeight;
    //@}

public:
    /** @name Constructors */
    //@{
    /** @brief Constructor with inputs
        @param r drum radius
        @param h drum height */
    __HOSTDEVICE__
    Drum(T r = T(0), T h = T(0));

    /** @brief Constructor with an input stream
        @param fileIn input stream */
    __HOST__
    Drum(std::istream& fileIn);

    /** @brief Constructor with an XML node as an input parameter
        @param root XML node */
    __HOST__
    Drum(DOMNode* root);

    /** @brief Destructor */
    __HOSTDEVICE__
    ~Drum();
    //@}

    /** @name Get methods */
    //@{
    /** @brief Gets the convex type */
    __HOSTDEVICE__
    ConvexType getConvexType() const final;

    /** @brief Gets the drum radius */
    __HOSTDEVICE__
    T getRadius() const;

    /** @brief Gets the drum height */
    __HOSTDEVICE__
    T getHeight() const;
    //@}

    /** @name Methods */
    //@{
    /** @brief Returns a clone of the drum */
    __HOSTDEVICE__
    Convex<T>* clone() const final;

    /** @brief Returns the drum shell area used as geometric measure */
    __HOSTDEVICE__
    T computeVolume() const final;

    /** @brief Computes the thin-shell inertia tensor
        @param inertia inertia tensor */
    __HOSTDEVICE__
    void computeInertia(T (&inertia)[3]) const final;

    /** @brief Returns the circumscribed radius of the drum */
    __HOSTDEVICE__
    T computeCircumscribedRadius() const final;

    /** @brief Returns the half-length of the bounding box fitted to the drum without
        considering the transformation */
    __HOSTDEVICE__
    Vector3<T> computeBoundingBox() const final;

    /** @brief Returns the tightest bounding cylinder fitted to the drum in body-local frame. */
    __HOSTDEVICE__
    Vector3<T> computeBoundingCylinder() const final;

    /** @brief Proxy support function of the matching finite cylinder.
        Analytical inside-drum contact is handled separately in CollisionDetection.
        @param v direction */
    __HOSTDEVICE__
    Vector3<T> support(const Vector3<T>& v) const final;

    /** @brief Eroded proxy support of the matching finite cylinder */
    __HOSTDEVICE__
    Vector3<T> support(const Vector3<T>& v, T crust) const final;
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

    /** @brief Returns the number of points to write the drum in a Paraview format. */
    __HOST__
    int numberOfPoints_PARAVIEW() const final;

    /** @brief Returns the number of elementary polytopes to write the drum in a Paraview format.
     */
    __HOST__
    int numberOfCells_PARAVIEW() const final;

    /** @brief Returns a list of points describing the drum in a Paraview format.
        @param transform geometric transformation
        @param translation additional center of mass translation */
    __HOST__
    std::list<Vector3<T>> writePoints_PARAVIEW(const Transform3<T>& transform,
                                               Vector3<T> const*    translation) const final;

    /** @brief Writes the connectivity of the drum in a Paraview format
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

    /** @brief Returns drum shape parameters: [radius, halfHeight, 0, 0, 0] */
    __HOSTDEVICE__
    void getShapeParameters(T (&params)[5]) const final;

    /** @brief Returns "Drum" */
    __HOST__
    std::string getConvexName() const final;

    /** @brief Writes the drum as an OBJ mesh (see Convex::writeOBJ) */
    __HOST__
    void writeOBJ(std::ostream& f, size_t& firstpoint_number) const final;
    //@}
};

#endif