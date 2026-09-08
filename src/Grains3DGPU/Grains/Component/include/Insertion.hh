#ifndef _INSERTION_HH_
#define _INSERTION_HH_

#include <variant>

#include "GrainsMemBuffer.hh"
#include "InsertionWindow.hh"
#include "Kinematics.hh"
#include "Quaternion.hh"
#include "ReaderXML.hh"
#include "RigidBody.hh"

/** @name Enumerations */
//@{
/** @brief Various supported insertion types */
enum InsertionType
{
    /** @brief Random insertion, i.e., randomly drawn values */
    RANDOMINSERTION,
    /** @brief File insertion, i.e., reading from a given file */
    FILEINSERTION,
    /** @brief Constant insertion, i.e., all values are set to a constant */
    CONSTANTINSERTION,
    /** @brief Array insertion, i.e., placed on a regular 3D grid */
    ARRAYINSERTION,
    /** @brief Default insertion, i.e., all values are zero */
    DEFAULTINSERTION
};

/** @brief Various supported insertion modes */
enum InsertionMode
{
    /** @brief Insertion at the initial time */
    INITIAL,
    /** @brief Insertion over the simulation time */
    OVERTIME
};

/** @brief Layout used for array insertion coordinates */
enum ArrayInsertionLayout
{
    /** @brief Cartesian grid generated from DX/DY/DZ and NX/NY/NZ */
    ARRAY_LAYOUT_RECTANGULAR,
    /** @brief Cylindrical grid generated from DR/DT/DH and NR/NT/NH */
    ARRAY_LAYOUT_CYLINDRICAL
};
//@}

/** @name Structs */
//@{
/** @brief Struct to hold the parameters for array insertion */
template <typename T>
struct ArrayInsertionData
{
    /** @brief Starting corner of the 3D grid */
    Vector3<T> origin;
    /** @brief Grid spacing in X, Y, Z */
    Vector3<T> spacing;
    /** @brief Local insertion frame X direction in world coordinates */
    Vector3<T> basisX;
    /** @brief Local insertion frame Y direction in world coordinates */
    Vector3<T> basisY;
    /** @brief Local insertion frame Z direction in world coordinates */
    Vector3<T> basisZ;
    /** @brief Number of grid points in X */
    uint nx;
    /** @brief Number of grid points in Y */
    uint ny;
    /** @brief Number of grid points in Z */
    uint nz;
    /** @brief Coordinate layout for this array insertion */
    ArrayInsertionLayout layout;
    /** @brief Total number of generated insertion slots */
    uint totalSlots;
    /** @brief Flat index of the next grid point to return */
    uint idx;
    /** @brief Radial (cylindrical) or linear (rectangular) growth per height layer.
        When positive, each layer gains floor(ih*DA/DR) extra rings or floor(iz*DA/DX|DY)
        extra columns/rows. Zero means standard layout with no expansion. */
    T slopeDA;
};

/** @brief info required for comping up with an insertion position. It can be either a value (T)
    that is used as the seed for random generator algorithm, a string (std::string) that is used as
    the pathToFile, a 3d vector (vector3<T>) for constant values, an ArrayInsertionData<T> for
    regular 3D grid insertion, and a value (0) in case the default insertion option is desired. */
template <typename T>
using InsertionInfo = std::
    variant<std::vector<InsertionWindow<T>>, std::ifstream, Vector3<T>, ArrayInsertionData<T>>;
//@}

// =================================================================================================
/** @brief The class Insertion.

    This class provides funtionalities to insert components in the simulation.

    @author A.YAZDANI - 2024 - Construction */
// =================================================================================================
template <typename T>
class Insertion
{
protected:
    /** @name Parameters */
    //@{
    /** \brief insertion type for position */
    InsertionType m_positionType;
    /** \brief insertion type for orientation */
    InsertionType m_orientationType;
    /** \brief insertion type for translational velocity */
    InsertionType m_translationalVelType;
    /** \brief insertion type for angular velocity */
    InsertionType m_angularVelType;
    /** \brief info required for coming up with an insertion position. */
    InsertionInfo<T> m_positionInsertionInfo;
    /** \brief info required for coming up with an insertion orientation.*/
    InsertionInfo<T> m_orientationInsertionInfo;
    /** \brief info required for coming up with an insertion velocity. */
    InsertionInfo<T> m_translationalVelInsertionInfo;
    /** \brief info required for coming up with an insertion omega. */
    InsertionInfo<T> m_angularVelInsertionInfo;
    /** \brief If insertion should be forced (true) or not (false) */
    bool m_forceInsertion;
    //@}

public:
    /**@name Contructors */
    //@{
    /** @brief Default constructor */
    __HOST__
    Insertion();

    /** @brief Constructor with XML element
        @param dn XML element */
    __HOST__
    Insertion(DOMNode* dn);

    /** @brief Destructor */
    __HOST__
    ~Insertion();
    //@}

    /** @name Get methods */
    //@{
    // TODO: INSERTION MODE AND TYPE?
    //@}

    /** @name Set methods */
    //@{
    /** @brief Set position insertion info (type determined by variant content)
        @param info Insertion info variant (windows, ifstream, or Vector3) */
    __HOST__
    void setPositionInsertionInfo(InsertionInfo<T>&& info);

    /** @brief Set orientation insertion info (type determined by variant content)
        @param info Insertion info variant (windows, ifstream, or Vector3) */
    __HOST__
    void setOrientationInsertionInfo(InsertionInfo<T>&& info);

    /** @brief Set translational velocity insertion info (type determined by variant content)
        @param info Insertion info variant (windows, ifstream, or Vector3) */
    __HOST__
    void setTranslationalVelInsertionInfo(InsertionInfo<T>&& info);

    /** @brief Set angular velocity insertion info (type determined by variant content)
        @param info Insertion info variant (windows, ifstream, or Vector3) */
    __HOST__
    void setAngularVelInsertionInfo(InsertionInfo<T>&& info);

    /** @brief Set force insertion flag
        @param forceInsertion If true, skip overlap checking during insertion */
    __HOST__
    void setForceInsertion(bool forceInsertion);
    //@}

    /** @name Methods */
    //@{
    /** @brief Reads an XML node to set the insertion type and related info
        for each component
        @param dn XML node */
    __HOST__
    std::pair<InsertionType, InsertionInfo<T>> readTypeAndData(DOMNode* root);

    /** @brief Returns a vector of Vector3 accroding to type and data.
        @param type insertion type
        @param data insertion info */
    __HOST__
    Vector3<T> fetchInsertionData(InsertionType const type, InsertionInfo<T>& data);

    /** @brief Returns all required data members to insert components as a vector.
        @param rigidBody rigid body buffer
        @param position position buffer
        @param quaternion quaternion buffer
        @param velocity velocity buffer
        @param LCParameters linked cell parameters
        @param numObstacles number of obstacles
        @param numParticles number of particles */
    __HOST__
    void insert(const GrainsMemBuffer<RigidBody<T>*>* rigidBody,
                GrainsMemBuffer<Vector3<T>>&          position,
                GrainsMemBuffer<Quaternion<T>>&       quaternion,
                GrainsMemBuffer<Kinematics<T>>&       velocity,
                const LinkedCellParameters<T>&        LCParameters,
                const uint                            numObstacles,
                const uint                            numParticles,
                const GrainsMemBuffer<uint>&          bodyTag,
                const GrainsMemBuffer<Vector3<T>>&    localPos,
                const GrainsMemBuffer<Quaternion<T>>& localQuat);
    //@}
};

#endif
