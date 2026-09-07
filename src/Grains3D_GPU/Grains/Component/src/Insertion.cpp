#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>

#include "BodyTag.hh"
#include "GJK.hh"
#include "GrainsMemBuffer.hh"
#include "GrainsUtils.hh"
#include "Insertion.hh"
#include "LinkedCell_Host.hh"
#include "OBB.hh"
#include "QuaternionMath.hh"

/* ============================================================================================== */
/* Low-Level Methods                                                                              */
/* ============================================================================================== */
// Helper function to build an insertion frame (basisX, basisY, basisZ) from a given Z direction
// vector
template <typename T>
__HOST__ static INLINE void buildInsertionFrameFromZ(const Vector3<T>& zDir,
                                                     Vector3<T>&       basisX,
                                                     Vector3<T>&       basisY,
                                                     Vector3<T>&       basisZ)
{
    GAssert(norm2(zDir) > T(0), "Insertion direction vector must be non-zero! Aborting Grains!");
    basisZ = zDir.normalized();

    // Keep the default rectangular frame unchanged when ZDirection is (0,0,1).
    Vector3<T> ref(T(0), T(1), T(0));
    if(std::abs(basisZ * ref) > T(0.95))
        ref = Vector3<T>(T(1), T(0), T(0));

    basisX = (ref ^ basisZ).normalized();
    basisY = basisZ ^ basisX;
}

// -------------------------------------------------------------------------------------------------
// Reads if the root is of type Random
template <typename T>
__HOST__ static INLINE InsertionInfo<T> readDataRand(DOMNode* root)
{
    // Random generator seed. We pass it to InsertionWindow directly.
    // We also set the seed with srand. We use it for to randomly pick an
    // insertion window
    RandomGeneratorSeed rgs;
    unsigned            cSeed      = 1;  // default C RNG seed for window selection
    std::string         seedString = ReaderXML::getNodeAttr_String(root, "Seed");
    if(seedString == "UserDefined")
    {
        uint val = ReaderXML::getNodeAttr_Int(root, "Value");
        GAssert(val, "Seed value is not provided. Aborting Grains!");
        rgs   = RGS_UDEF;
        cSeed = static_cast<unsigned>(val);
        GoutWI(12, "Random initialization with seed", std::to_string(val) + ".");
    }
    else if(seedString == "Random")
    {
        rgs   = RGS_RANDOM;
        cSeed = static_cast<unsigned>(time(NULL));
        GoutWI(12, "Random initialization with random seed.");
    }
    // if ( seedString == "Default" )
    else
    {
        rgs   = RGS_DEFAULT;
        cSeed = 1u;
        GoutWI(12, "Random initialization with default seed.");
    }
    // Seed C RNG for selecting among multiple insertion windows
    srand(cSeed);

    // Insertion window
    DOMNode*                        nWindows = ReaderXML::getNode(root, "Windows");
    std::vector<InsertionWindow<T>> insertionWindows;
    if(nWindows)
    {
        // DOMNodeList* allWindows = ReaderXML::getNodes( nWindows, "Window" );
        DOMNodeList* allWindows = ReaderXML::getNodes(nWindows);
        for(int i = 0; i < allWindows->getLength(); i++)
        {
            DOMNode* nWindow = allWindows->item(i);
            insertionWindows.push_back(InsertionWindow<T>(nWindow, rgs));
        }
    }

    return (insertionWindows);
}

// -------------------------------------------------------------------------------------------------
// Reads if the root is of type File
template <typename T>
__HOST__ static INLINE InsertionInfo<T> readDataFile(DOMNode* root)
{
    std::string   fileName = ReaderXML::getNodeAttr_String(root, "Name");
    std::ifstream file(fileName);
    GAssert(file.good(), "File initialization failed! Aborting Grains!");
    GoutWI(12, "File initialization with path" + fileName + ".");

    return (file);
}

// -------------------------------------------------------------------------------------------------
// Reads if the root is of type Constant
template <typename T>
__HOST__ static INLINE InsertionInfo<T> readDataCons(DOMNode* root)
{
    T          xVal = T(ReaderXML::getNodeAttr_Double(root, "X"));
    T          yVal = T(ReaderXML::getNodeAttr_Double(root, "Y"));
    T          zVal = T(ReaderXML::getNodeAttr_Double(root, "Z"));
    Vector3<T> vec(xVal, yVal, zVal);
    GoutWI(12, "Constant initialization with", Vector3ToString(vec), ".");

    return (vec);
}

// -------------------------------------------------------------------------------------------------
// Reads if the root is of type Zero
template <typename T>
__HOST__ static INLINE InsertionInfo<T> readDataZero(DOMNode* root)
{
    Vector3<T> vec(T(0), T(0), T(0));
    GoutWI(12, "Zero initialization.");

    return (vec);
}

// -------------------------------------------------------------------------------------------------
// Reads if the root is of type Array (regular 3D grid)
template <typename T>
__HOST__ static INLINE InsertionInfo<T> readDataArray(DOMNode* root)
{
    ArrayInsertionData<T> aid;

    DOMNode* nOrigin  = ReaderXML::getNode(root, "Origin");
    DOMNode* nSpacing = ReaderXML::getNode(root, "Spacing");
    DOMNode* nCount   = ReaderXML::getNode(root, "Count");
    GAssert(nOrigin, "Array insertion requires an <Origin> child node! Aborting Grains!");
    GAssert(nSpacing, "Array insertion requires a <Spacing> child node! Aborting Grains!");
    GAssert(nCount, "Array insertion requires a <Count> child node! Aborting Grains!");

    aid.origin.setValue(T(ReaderXML::getNodeAttr_Double(nOrigin, "X")),
                        T(ReaderXML::getNodeAttr_Double(nOrigin, "Y")),
                        T(ReaderXML::getNodeAttr_Double(nOrigin, "Z")));

    const bool hasDX = ReaderXML::hasNodeAttr(nSpacing, "DX");
    const bool hasDY = ReaderXML::hasNodeAttr(nSpacing, "DY");
    const bool hasDZ = ReaderXML::hasNodeAttr(nSpacing, "DZ");
    const bool hasDR = ReaderXML::hasNodeAttr(nSpacing, "DR");
    const bool hasDT = ReaderXML::hasNodeAttr(nSpacing, "DT");
    const bool hasDH = ReaderXML::hasNodeAttr(nSpacing, "DH");

    const bool hasNX = ReaderXML::hasNodeAttr(nCount, "NX");
    const bool hasNY = ReaderXML::hasNodeAttr(nCount, "NY");
    const bool hasNZ = ReaderXML::hasNodeAttr(nCount, "NZ");
    const bool hasNR = ReaderXML::hasNodeAttr(nCount, "NR");
    const bool hasNT = ReaderXML::hasNodeAttr(nCount, "NT");
    const bool hasNH = ReaderXML::hasNodeAttr(nCount, "NH");

    const bool rectSpacingComplete = hasDX && hasDY && hasDZ;
    const bool cylSpacingComplete  = hasDR && hasDT && hasDH;
    const bool rectCountComplete   = hasNX && hasNY && hasNZ;
    const bool cylCountComplete    = hasNR && hasNT && hasNH;

    GAssert(rectSpacingComplete || cylSpacingComplete,
            "Array insertion <Spacing> must provide either DX/DY/DZ (rectangular) or "
            "DR/DT/DH (cylindrical)! Aborting Grains!");
    GAssert(!(rectSpacingComplete && cylSpacingComplete),
            "Array insertion <Spacing> cannot mix rectangular and cylindrical attributes! "
            "Aborting Grains!");
    GAssert(rectCountComplete || cylCountComplete,
            "Array insertion <Count> must provide either NX/NY/NZ (rectangular) or "
            "NR/NT/NH (cylindrical)! Aborting Grains!");
    GAssert(!(rectCountComplete && cylCountComplete),
            "Array insertion <Count> cannot mix rectangular and cylindrical attributes! "
            "Aborting Grains!");

    GAssert(rectSpacingComplete == rectCountComplete,
            "Array insertion mode mismatch: DX/DY/DZ must be paired with NX/NY/NZ, and "
            "DR/DT/DH must be paired with NR/NT/NH! Aborting Grains!");

    DOMNode* nZDirection = ReaderXML::getNode(root, "ZDirection");

    if(rectSpacingComplete)
    {
        aid.layout = ARRAY_LAYOUT_RECTANGULAR;
        aid.spacing.setValue(T(ReaderXML::getNodeAttr_Double(nSpacing, "DX")),
                             T(ReaderXML::getNodeAttr_Double(nSpacing, "DY")),
                             T(ReaderXML::getNodeAttr_Double(nSpacing, "DZ")));
        aid.nx      = static_cast<uint>(ReaderXML::getNodeAttr_Int(nCount, "NX"));
        aid.ny      = static_cast<uint>(ReaderXML::getNodeAttr_Int(nCount, "NY"));
        aid.nz      = static_cast<uint>(ReaderXML::getNodeAttr_Int(nCount, "NZ"));
        aid.slopeDA = ReaderXML::hasNodeAttr(nSpacing, "DA")
                          ? T(ReaderXML::getNodeAttr_Double(nSpacing, "DA"))
                          : T(0);
        if(aid.slopeDA > T(0))
        {
            aid.totalSlots = 0u;
            for(uint iz = 0u; iz < aid.nz; ++iz)
            {
                const uint nxLayer
                    = aid.nx + static_cast<uint>(T(iz) * aid.slopeDA / aid.spacing[0]);
                const uint nyLayer
                    = aid.ny + static_cast<uint>(T(iz) * aid.slopeDA / aid.spacing[1]);
                aid.totalSlots += nxLayer * nyLayer;
            }
        }
        else
            aid.totalSlots = aid.nx * aid.ny * aid.nz;

        Vector3<T> zDir(T(0), T(0), T(1));
        if(nZDirection)
            zDir.setValue(T(ReaderXML::getNodeAttr_Double(nZDirection, "X")),
                          T(ReaderXML::getNodeAttr_Double(nZDirection, "Y")),
                          T(ReaderXML::getNodeAttr_Double(nZDirection, "Z")));
        buildInsertionFrameFromZ(zDir, aid.basisX, aid.basisY, aid.basisZ);
    }
    else
    {
        aid.spacing.setValue(T(ReaderXML::getNodeAttr_Double(nSpacing, "DR")),
                             T(ReaderXML::getNodeAttr_Double(nSpacing, "DT")),
                             T(ReaderXML::getNodeAttr_Double(nSpacing, "DH")));
        aid.nx = static_cast<uint>(ReaderXML::getNodeAttr_Int(nCount, "NR"));
        aid.ny = static_cast<uint>(ReaderXML::getNodeAttr_Int(nCount, "NT"));
        aid.nz = static_cast<uint>(ReaderXML::getNodeAttr_Int(nCount, "NH"));

        GAssert(aid.spacing[0] > T(0) && aid.spacing[1] > T(0) && aid.spacing[2] > T(0),
                "Cylindrical array insertion requires positive DR/DT/DH! Aborting Grains!");

        aid.layout  = ARRAY_LAYOUT_CYLINDRICAL;
        aid.slopeDA = ReaderXML::hasNodeAttr(nSpacing, "DA")
                          ? T(ReaderXML::getNodeAttr_Double(nSpacing, "DA"))
                          : T(0);

        const T twoPi  = T(2) * T(std::acos(-1.0));
        aid.totalSlots = 0u;
        for(uint ih = 0u; ih < aid.nz; ++ih)
        {
            // When DA is provided, each height layer gains floor(ih*DA/DR) extra rings.
            const uint nrLayer = aid.nx + static_cast<uint>(T(ih) * aid.slopeDA / aid.spacing[0]);
            for(uint ir = 0u; ir < nrLayer; ++ir)
            {
                uint ringNt = 1u;
                if(ir > 0u)
                {
                    const T r         = T(ir) * aid.spacing[0];
                    const T minDTheta = aid.spacing[0] / r;
                    const T dThetaEff = std::max(aid.spacing[1], minDTheta);
                    ringNt            = std::max(1u, static_cast<uint>(twoPi / dThetaEff));
                    ringNt            = std::min(aid.ny, ringNt);
                }
                aid.totalSlots += ringNt;
            }
        }

        Vector3<T> zDir(T(0), T(1), T(0));
        if(nZDirection)
            zDir.setValue(T(ReaderXML::getNodeAttr_Double(nZDirection, "X")),
                          T(ReaderXML::getNodeAttr_Double(nZDirection, "Y")),
                          T(ReaderXML::getNodeAttr_Double(nZDirection, "Z")));
        buildInsertionFrameFromZ(zDir, aid.basisX, aid.basisY, aid.basisZ);
    }
    aid.idx = 0;

    GAssert(aid.nx > 0 && aid.ny > 0 && aid.nz > 0,
            "Array insertion counts must all be positive! Aborting Grains!");

    if(aid.layout == ARRAY_LAYOUT_RECTANGULAR)
    {
        if(aid.slopeDA > T(0))
            GoutWI(12,
                   "Array initialization (rectangular, DA=" + std::to_string(aid.slopeDA) + "):",
                   std::to_string(aid.nx) + "x" + std::to_string(aid.ny) + "x"
                       + std::to_string(aid.nz),
                   "base,",
                   std::to_string(aid.totalSlots),
                   "effective grid points.");
        else
            GoutWI(12,
                   "Array initialization (rectangular):",
                   std::to_string(aid.nx) + "x" + std::to_string(aid.ny) + "x"
                       + std::to_string(aid.nz),
                   "=",
                   std::to_string(aid.totalSlots),
                   "grid points.");
    }
    else
    {
        if(aid.slopeDA > T(0))
            GoutWI(12,
                   "Array initialization (cylindrical, adaptive theta, DA="
                       + std::to_string(aid.slopeDA) + "):",
                   std::to_string(aid.nx) + "x" + std::to_string(aid.ny) + "x"
                       + std::to_string(aid.nz),
                   "base,",
                   std::to_string(aid.totalSlots),
                   "effective grid points.");
        else
            GoutWI(12,
                   "Array initialization (cylindrical, adaptive theta):",
                   std::to_string(aid.nx) + "x" + std::to_string(aid.ny) + "x"
                       + std::to_string(aid.nz),
                   "requested,",
                   std::to_string(aid.totalSlots),
                   "effective grid points.");
    }

    return (aid);
}

/* ============================================================================================== */
/* High-Level Methods                                                                             */
/* ============================================================================================== */
// Default constructor
template <typename T>
__HOST__ Insertion<T>::Insertion()
    : m_positionType(DEFAULTINSERTION)
    , m_orientationType(DEFAULTINSERTION)
    , m_translationalVelType(DEFAULTINSERTION)
    , m_angularVelType(DEFAULTINSERTION)
    , m_positionInsertionInfo(Vector3<T>())
    , m_orientationInsertionInfo(Vector3<T>())
    , m_translationalVelInsertionInfo(Vector3<T>())
    , m_angularVelInsertionInfo(Vector3<T>())
    , m_forceInsertion(false)
{
}

// -------------------------------------------------------------------------------------------------
// Constructor with XML node
template <typename T>
__HOST__ Insertion<T>::Insertion(DOMNode* dn)
{
    // We define a lambda function to read the XML node
    auto read = [](DOMNode* root, InsertionType& type, InsertionInfo<T>& data) {
        std::string nType = ReaderXML::getNodeAttr_String(root, "Type");
        if(nType == "Random")
        {
            type = RANDOMINSERTION;
            data = readDataRand<T>(root);
        }
        else if(nType == "File")
        {
            type = FILEINSERTION;
            data = readDataFile<T>(root);
        }
        else if(nType == "Constant")
        {
            type = CONSTANTINSERTION;
            data = readDataCons<T>(root);
        }
        else if(nType == "Array")
        {
            type = ARRAYINSERTION;
            data = readDataArray<T>(root);
        }
        else if(nType == "Zero")
        {
            type = DEFAULTINSERTION;
            data = readDataZero<T>(root);
        }
        else
            GAbort("Unknown Type in ParticleInsertion! Aborting Grains!");
    };

    GAssert(dn, "ParticleInsertion node is missing! Aborting Grains!");

    GoutWI(9, "Reading PositionInsertion Policy ...");
    if(ReaderXML::getNode(dn, "InitialPosition"))
    {
        DOMNode* nIP = ReaderXML::getNode(dn, "InitialPosition");
        read(nIP, m_positionType, m_positionInsertionInfo);
    }
    else
    {
        GAbort("InitialPosition node is missing in ParticleInsertion! Aborting "
               "Grains!");
    }

    GoutWI(9, "Reading OrientationInsertion Policy ...");
    if(ReaderXML::getNode(dn, "InitialOrientation"))
    {
        DOMNode* nIO = ReaderXML::getNode(dn, "InitialOrientation");
        read(nIO, m_orientationType, m_orientationInsertionInfo);
    }
    else
    {
        m_orientationType          = DEFAULTINSERTION;
        m_orientationInsertionInfo = Vector3<T>(T(0), T(0), T(0));
        GoutWI(12, "No InitialOrientation node found. Using default.");
    }

    GoutWI(9, "Reading VeclocityInsertion Policy ...");
    if(ReaderXML::getNode(dn, "InitialVelocity"))
    {
        DOMNode* nIV = ReaderXML::getNode(dn, "InitialVelocity");
        read(nIV, m_translationalVelType, m_translationalVelInsertionInfo);
    }
    else
    {
        m_translationalVelType          = DEFAULTINSERTION;
        m_translationalVelInsertionInfo = Vector3<T>(T(0), T(0), T(0));
        GoutWI(12, "No InitialVelocity node found. Using default.");
    }

    GoutWI(9, "Reading AngularVeclocityInsertion Policy ...");
    if(ReaderXML::getNode(dn, "InitialAngularVelocity"))
    {
        DOMNode* nIA = ReaderXML::getNode(dn, "InitialAngularVelocity");
        read(nIA, m_angularVelType, m_angularVelInsertionInfo);
    }
    else
    {
        m_angularVelType          = DEFAULTINSERTION;
        m_angularVelInsertionInfo = Vector3<T>(T(0), T(0), T(0));
        GoutWI(12, "No InitialAngularVelocity node found. Using default.");
    }

    if(ReaderXML::hasNodeAttr(dn, "ForceInsertion"))
        m_forceInsertion = static_cast<bool>(ReaderXML::getNodeAttr_Int(dn, "ForceInsertion"));
    else
        m_forceInsertion = false;
}

// -------------------------------------------------------------------------------------------------
// Destructor
template <typename T>
__HOST__ Insertion<T>::~Insertion()
{
    if(std::holds_alternative<std::ifstream>(m_positionInsertionInfo))
        (std::get<std::ifstream>(m_positionInsertionInfo)).close();
    if(std::holds_alternative<std::ifstream>(m_orientationInsertionInfo))
        (std::get<std::ifstream>(m_orientationInsertionInfo)).close();
    if(std::holds_alternative<std::ifstream>(m_translationalVelInsertionInfo))
        (std::get<std::ifstream>(m_translationalVelInsertionInfo)).close();
    if(std::holds_alternative<std::ifstream>(m_angularVelInsertionInfo))
        (std::get<std::ifstream>(m_angularVelInsertionInfo)).close();
}

// -------------------------------------------------------------------------------------------------
// Set position insertion info (type determined by variant content)
template <typename T>
__HOST__ void Insertion<T>::setPositionInsertionInfo(InsertionInfo<T>&& info)
{
    // Determine type from variant content BEFORE moving
    if(std::holds_alternative<std::vector<InsertionWindow<T>>>(info))
        m_positionType = RANDOMINSERTION;
    else if(std::holds_alternative<std::ifstream>(info))
        m_positionType = FILEINSERTION;
    else if(std::holds_alternative<ArrayInsertionData<T>>(info))
        m_positionType = ARRAYINSERTION;
    else if(std::holds_alternative<Vector3<T>>(info))
        m_positionType = CONSTANTINSERTION;
    else
        m_positionType = DEFAULTINSERTION;

    m_positionInsertionInfo = std::move(info);
}

// -------------------------------------------------------------------------------------------------
// Set orientation insertion info (type determined by variant content)
template <typename T>
__HOST__ void Insertion<T>::setOrientationInsertionInfo(InsertionInfo<T>&& info)
{
    if(std::holds_alternative<std::vector<InsertionWindow<T>>>(info))
        m_orientationType = RANDOMINSERTION;
    else if(std::holds_alternative<std::ifstream>(info))
        m_orientationType = FILEINSERTION;
    else if(std::holds_alternative<ArrayInsertionData<T>>(info))
        m_orientationType = ARRAYINSERTION;
    else if(std::holds_alternative<Vector3<T>>(info))
        m_orientationType = CONSTANTINSERTION;
    else
        m_orientationType = DEFAULTINSERTION;

    m_orientationInsertionInfo = std::move(info);
}

// -------------------------------------------------------------------------------------------------
// Set translational velocity insertion info (type determined by variant content)
template <typename T>
__HOST__ void Insertion<T>::setTranslationalVelInsertionInfo(InsertionInfo<T>&& info)
{
    if(std::holds_alternative<std::vector<InsertionWindow<T>>>(info))
        m_translationalVelType = RANDOMINSERTION;
    else if(std::holds_alternative<std::ifstream>(info))
        m_translationalVelType = FILEINSERTION;
    else if(std::holds_alternative<ArrayInsertionData<T>>(info))
        m_translationalVelType = ARRAYINSERTION;
    else if(std::holds_alternative<Vector3<T>>(info))
        m_translationalVelType = CONSTANTINSERTION;
    else
        m_translationalVelType = DEFAULTINSERTION;

    m_translationalVelInsertionInfo = std::move(info);
}

// -------------------------------------------------------------------------------------------------
// Set angular velocity insertion info (type determined by variant content)
template <typename T>
__HOST__ void Insertion<T>::setAngularVelInsertionInfo(InsertionInfo<T>&& info)
{
    if(std::holds_alternative<std::vector<InsertionWindow<T>>>(info))
        m_angularVelType = RANDOMINSERTION;
    else if(std::holds_alternative<std::ifstream>(info))
        m_angularVelType = FILEINSERTION;
    else if(std::holds_alternative<ArrayInsertionData<T>>(info))
        m_angularVelType = ARRAYINSERTION;
    else if(std::holds_alternative<Vector3<T>>(info))
        m_angularVelType = CONSTANTINSERTION;
    else
        m_angularVelType = DEFAULTINSERTION;

    m_angularVelInsertionInfo = std::move(info);
}

// -------------------------------------------------------------------------------------------------
// Set force insertion flag
template <typename T>
__HOST__ void Insertion<T>::setForceInsertion(bool forceInsertion)
{
    m_forceInsertion = forceInsertion;
}

// -------------------------------------------------------------------------------------------------
// Returns a vector of Vector3 accroding to type and data
template <typename T>
__HOST__ Vector3<T> Insertion<T>::fetchInsertionData(InsertionType const type,
                                                     InsertionInfo<T>&   data)
{
    // We only return a vector3. It is clear how it works for position, and
    // kinematics. However, for orientation, it returns the vector3 of rotation
    // angles. We later construct a quaternion.
    if(type == RANDOMINSERTION)
    {
        auto& IWs = std::get<std::vector<InsertionWindow<T>>>(data);
        GAssert(!IWs.empty(), "Random insertion selected but no InsertionWindow defined!");
        if(IWs.size() == 1)
            return IWs[0].generateRandomPoint();

        // Randomly choose between the available insertion windows
        int random_IW = static_cast<int>(rand() % IWs.size());
        return IWs[random_IW].generateRandomPoint();
    }
    else if(type == FILEINSERTION)
    {
        Vector3<T> output;
        std::get<std::ifstream>(data) >> output;
        return (output);
    }
    else if(type == ARRAYINSERTION)
    {
        auto& aid = std::get<ArrayInsertionData<T>>(data);
        GAssert(aid.idx < aid.totalSlots,
                "Array insertion exhausted all grid slots! Increase counts in <Count>.");
        const uint flat = aid.idx++;

        // Generate small random perturbations (~1e-10)
        const T epsilon = T(1e-10);
        const T randX   = epsilon * T(2.0 * (double)rand() / RAND_MAX - 1.0);
        const T randY   = epsilon * T(2.0 * (double)rand() / RAND_MAX - 1.0);
        const T randZ   = epsilon * T(2.0 * (double)rand() / RAND_MAX - 1.0);

        if(aid.layout == ARRAY_LAYOUT_CYLINDRICAL)
        {
            // Cylindrical coordinates with adaptive angular sampling.
            // When slopeDA > 0 each height layer gains floor(ih*DA/DR) extra rings,
            // producing an expanding-cone packing.
            const T twoPi   = T(2) * T(std::acos(-1.0));
            uint    rem     = flat;
            uint    ih      = 0u;
            uint    ir      = 0u;
            uint    it      = 0u;
            uint    ihFound = 0u;

            bool found = false;
            for(ih = 0u; ih < aid.nz && !found; ++ih)
            {
                const uint nrLayer
                    = aid.nx + static_cast<uint>(T(ih) * aid.slopeDA / aid.spacing[0]);
                for(ir = 0u; ir < nrLayer; ++ir)
                {
                    uint ringNt = 1u;
                    if(ir > 0u)
                    {
                        const T r         = T(ir) * aid.spacing[0];
                        const T minDTheta = aid.spacing[0] / r;
                        const T dThetaEff = std::max(aid.spacing[1], minDTheta);
                        ringNt            = std::max(1u, static_cast<uint>(twoPi / dThetaEff));
                        ringNt            = std::min(aid.ny, ringNt);
                    }
                    if(rem < ringNt)
                    {
                        it      = rem;
                        ihFound = ih;
                        found   = true;
                        break;
                    }
                    rem -= ringNt;
                }
            }

            GAssert(found, "Internal error mapping cylindrical insertion index. Aborting Grains!");
            const T r = T(ir) * aid.spacing[0];
            T       theta(0);
            if(ir > 0u)
            {
                const T minDTheta = aid.spacing[0] / r;
                const T dThetaEff = std::max(aid.spacing[1], minDTheta);
                theta             = T(it) * dThetaEff;
            }
            const T h  = T(ihFound) * aid.spacing[2];
            const T lx = r * T(std::cos(theta));
            const T ly = r * T(std::sin(theta));
            const T lz = h;
            return Vector3<T>(aid.origin[0] + lx * aid.basisX[0] + ly * aid.basisY[0]
                                  + lz * aid.basisZ[0] + randX,
                              aid.origin[1] + lx * aid.basisX[1] + ly * aid.basisY[1]
                                  + lz * aid.basisZ[1] + randY,
                              aid.origin[2] + lx * aid.basisX[2] + ly * aid.basisY[2]
                                  + lz * aid.basisZ[2] + randZ);
        }

        if(aid.slopeDA > T(0))
        {
            // Rectangular grid expanding outward by DA per DZ layer.
            uint rem     = flat;
            uint ixFound = 0u;
            uint iyFound = 0u;
            uint izFound = 0u;
            bool found   = false;
            for(uint iz = 0u; iz < aid.nz && !found; ++iz)
            {
                const uint nxLayer
                    = aid.nx + static_cast<uint>(T(iz) * aid.slopeDA / aid.spacing[0]);
                const uint nyLayer
                    = aid.ny + static_cast<uint>(T(iz) * aid.slopeDA / aid.spacing[1]);
                const uint layerSlots = nxLayer * nyLayer;
                if(rem < layerSlots)
                {
                    ixFound = rem % nxLayer;
                    iyFound = rem / nxLayer;
                    izFound = iz;
                    found   = true;
                }
                else
                    rem -= layerSlots;
            }
            GAssert(found, "Internal error mapping rectangular insertion index. Aborting Grains!");
            const T lx = T(ixFound) * aid.spacing[0];
            const T ly = T(iyFound) * aid.spacing[1];
            const T lz = T(izFound) * aid.spacing[2];
            return Vector3<T>(aid.origin[0] + lx * aid.basisX[0] + ly * aid.basisY[0]
                                  + lz * aid.basisZ[0] + randX,
                              aid.origin[1] + lx * aid.basisX[1] + ly * aid.basisY[1]
                                  + lz * aid.basisZ[1] + randY,
                              aid.origin[2] + lx * aid.basisX[2] + ly * aid.basisY[2]
                                  + lz * aid.basisZ[2] + randZ);
        }

        const uint ix = flat % aid.nx;
        const uint iy = (flat / aid.nx) % aid.ny;
        const uint iz = flat / (aid.nx * aid.ny);

        const T lx = T(ix) * aid.spacing[0];
        const T ly = T(iy) * aid.spacing[1];
        const T lz = T(iz) * aid.spacing[2];
        return Vector3<T>(
            aid.origin[0] + lx * aid.basisX[0] + ly * aid.basisY[0] + lz * aid.basisZ[0] + randX,
            aid.origin[1] + lx * aid.basisX[1] + ly * aid.basisY[1] + lz * aid.basisZ[1] + randY,
            aid.origin[2] + lx * aid.basisX[2] + ly * aid.basisY[2] + lz * aid.basisZ[2] + randZ);
    }
    else if(type == CONSTANTINSERTION)
        return (std::get<Vector3<T>>(data));
    else
        return (Vector3<T>());
}

// -------------------------------------------------------------------------------------------------
// Populates position, orientation, and kinematics according to the insertion
// policy
template <typename T>
__HOST__ void Insertion<T>::insert(const GrainsMemBuffer<RigidBody<T>*>* rigidBody,
                                   GrainsMemBuffer<Vector3<T>>&          position,
                                   GrainsMemBuffer<Quaternion<T>>&       orientation,
                                   GrainsMemBuffer<Kinematics<T>>&       kinematics,
                                   const LinkedCellParameters<T>&        LCParameters,
                                   const uint                            numObstacles,
                                   const uint                            numParticles,
                                   const GrainsMemBuffer<uint>&          bodyTag,
                                   const GrainsMemBuffer<Vector3<T>>&    localPos,
                                   const GrainsMemBuffer<Quaternion<T>>& localQuat)
{
    Gout("Inserting", numParticles, "particles ...");

    // Progress bar: report every 5% of numParticles slots processed
    const uint progressInterval = std::max(1u, numParticles / 20u);
    uint       nextMilestone    = progressInterval;
    auto       reportProgress   = [&](uint i) {
        if(i + 1u >= nextMilestone || i + 1u == numParticles)
        {
            int pct = static_cast<int>(100u * (i + 1u) / numParticles);
            GoutWI(3,
                   "[Insertion]",
                   pct,
                   "% (" + std::to_string(i + 1u) + " / " + std::to_string(numParticles)
                       + " particles)");
            nextMilestone += progressInterval;
        }
    };

    if(m_forceInsertion || m_positionType == ARRAYINSERTION)
    {
        for(uint i = 0; i < numParticles; ++i)
        {
            const uint insertID = i + numObstacles;

            // Slave sub-bodies are positioned by updateSubBodyPositions() after the master is
            // placed -- skip them here to avoid consuming insertion data for non-master slots
            if(isSubBody(bodyTag[insertID]) && getSubBodyLocalIdx(bodyTag[insertID]) > 0u)
                continue;

            position[insertID] = fetchInsertionData(m_positionType, m_positionInsertionInfo);

            // Orientation angles. These are not matrices, so we have to
            // compute the quaternions later.
            Vector3<T> ori = fetchInsertionData(m_orientationType, m_orientationInsertionInfo);
            orientation[insertID] = Quaternion<T>(ori[X], ori[Y], ori[Z]) * orientation[insertID];

            Vector3<T> vel
                = fetchInsertionData(m_translationalVelType, m_translationalVelInsertionInfo);

            Vector3<T> ang       = fetchInsertionData(m_angularVelType, m_angularVelInsertionInfo);
            kinematics[insertID] = Kinematics<T>(vel, ang);
            reportProgress(i);
        }
    }
    else
    {
        // Max attempts to place a particle
        const uint maxAttempts = 1000;

        // Build a temporary linked-cell structure for strict insertion checks
        LinkedCell_Host<T> LC(rigidBody,
                              position,
                              orientation,
                              LCParameters,
                              numObstacles,
                              numParticles);

        // Helper: test one body against all already-placed particles via the LC.
        // Returns true if the body overlaps any existing particle.
        auto overlapsLC = [&](const uint           bodyID,
                              const Convex<T>&     convexTest,
                              const T              crustTest,
                              const Vector3<T>&    worldPos,
                              const Quaternion<T>& worldQuat) -> bool {
            std::vector<uint> neighborList;
            LC.collectPotentialNeighbors(worldPos, bodyID, neighborList);
            for(uint j : neighborList)
            {
                const RigidBody<T>* rbJ     = (*rigidBody)[j];
                const Convex<T>&    convexJ = *rbJ->getConvex();
                bool BVintersect = intersectOrientedBoundingBox(convexJ.computeBoundingBox(),
                                                                convexTest.computeBoundingBox(),
                                                                position[j],
                                                                worldPos,
                                                                orientation[j],
                                                                worldQuat);
                if(BVintersect)
                {
                    Vector3<T> pa, pb;
                    uint       nbIter = 0;
                    const T    crustJ = rbJ->getCrustThickness();
                    const T    gap    = computeClosestPoints_GJK<T, GJKType::JOHNSON>(convexJ,
                                                                                convexTest,
                                                                                position[j],
                                                                                worldPos,
                                                                                orientation[j],
                                                                                worldQuat,
                                                                                crustJ,
                                                                                crustTest,
                                                                                pa,
                                                                                pb,
                                                                                nbIter);
                    // computeClosestPoints_GJK returns G + crustA + crustB where G is the
                    // actual geometric gap.  Contact in the simulation is detected when
                    // G <= 0, i.e. gap <= crustA + crustB.
                    if(gap < crustJ + crustTest)
                        return true;
                }
            }
            return false;
        };

        // Overlap test: checks the master (or standalone) convex plus all slave sub-bodies
        // at their rigid-body offsets from the candidate master pose.
        auto canInsert = [&](const uint           insertID,
                             const Vector3<T>&    insertPosition,
                             const Quaternion<T>& insertQuaternion) -> bool {
            // Check the master / standalone body itself
            const RigidBody<T>* rbMaster     = (*rigidBody)[insertID];
            const Convex<T>&    convexMaster = *rbMaster->getConvex();
            const T             crustMaster  = rbMaster->getCrustThickness();

            // For a composite master (LocalIdx==0), insertPosition is the composite CM;
            // slave's world position is CM + q >> localPos[*].
            const uint       tag          = bodyTag[insertID];
            const bool       isCompMaster = isSubBody(tag) && (getSubBodyLocalIdx(tag) == 0u);
            const Vector3<T> masterWorldPos
                = isCompMaster ? insertPosition + (insertQuaternion >> localPos[insertID])
                               : insertPosition;
            const Quaternion<T> masterWorldQuat
                = isCompMaster ? insertQuaternion * localQuat[insertID] : insertQuaternion;

            if(overlapsLC(insertID, convexMaster, crustMaster, masterWorldPos, masterWorldQuat))
                return false;

            // For composite masters: check every slave sub-body at its world pose
            if(isCompMaster)
            {
                const uint cIdx = getCompositeIdx(tag);
                for(uint k = insertID + 1;
                    k < numObstacles + numParticles && isSubBody(bodyTag[k])
                    && getCompositeIdx(bodyTag[k]) == cIdx && getSubBodyLocalIdx(bodyTag[k]) > 0u;
                    ++k)
                {
                    const Vector3<T> slavePos = insertPosition + (insertQuaternion >> localPos[k]);
                    const Quaternion<T> slaveQuat = insertQuaternion * localQuat[k];
                    const RigidBody<T>* rbK       = (*rigidBody)[k];
                    const Convex<T>&    convexK   = *rbK->getConvex();
                    if(overlapsLC(k, convexK, rbK->getCrustThickness(), slavePos, slaveQuat))
                        return false;
                }
            }
            return true;
        };

        // Inserting particles
        for(uint i = 0; i < numParticles; ++i)
        {
            const uint insertID = i + numObstacles;

            // Slave sub-bodies are positioned after their master is placed; skip them here
            if(isSubBody(bodyTag[insertID]) && getSubBodyLocalIdx(bodyTag[insertID]) > 0u)
                continue;

            bool placed = false;
            for(uint attempt = 0; attempt < maxAttempts && !placed; ++attempt)
            {
                const Vector3<T>& pCand
                    = fetchInsertionData(m_positionType, m_positionInsertionInfo);
                // Orientation angles. These are not matrices, so we have to
                // compute the quaternions later.
                Vector3<T> ori = fetchInsertionData(m_orientationType, m_orientationInsertionInfo);
                Quaternion<T>        quat(ori[X], ori[Y], ori[Z]);
                const Quaternion<T>& qCand = quat * orientation[insertID];
                // Check if candidate position is within domain bounds
                // clang-format off
                bool withinBounds = (pCand[0] >= LCParameters.minCorner[0] && 
                                     pCand[0] <= LCParameters.maxCorner[0] && 
                                     pCand[1] >= LCParameters.minCorner[1] && 
                                     pCand[1] <= LCParameters.maxCorner[1] &&
                                     pCand[2] >= LCParameters.minCorner[2] &&
                                     pCand[2] <= LCParameters.maxCorner[2]);
                // clang-format on
                if(withinBounds && canInsert(insertID, pCand, qCand))
                {
                    // pCand is the composite CM (or standalone body position).
                    // Store it directly so initCompositeFrames() can later snapshot the CM.
                    position[insertID]    = pCand;
                    orientation[insertID] = qCand;
                    Vector3<T> vel        = fetchInsertionData(m_translationalVelType,
                                                        m_translationalVelInsertionInfo);
                    Vector3<T> ang
                        = fetchInsertionData(m_angularVelType, m_angularVelInsertionInfo);
                    kinematics[insertID] = Kinematics<T>(vel, ang);
                    placed               = true;
                    // Add master to LC using its actual world position (CM + q >> localPos[0])
                    const Cells<T>*  cells   = LC.getLinkedCell()[0];
                    const uint       tag     = bodyTag[insertID];
                    const bool       isCompM = isSubBody(tag) && (getSubBodyLocalIdx(tag) == 0u);
                    const Vector3<T> masterWorldPos
                        = isCompM ? pCand + (qCand >> localPos[insertID]) : pCand;
                    LC.addParticleToCell(insertID, cells->computeCellHash(masterWorldPos));
                    // For composite masters: write all slave world poses and add to LC
                    if(isCompM)
                    {
                        const uint cIdx = getCompositeIdx(tag);
                        for(uint k = insertID + 1;
                            k < numObstacles + numParticles && isSubBody(bodyTag[k])
                            && getCompositeIdx(bodyTag[k]) == cIdx
                            && getSubBodyLocalIdx(bodyTag[k]) > 0u;
                            ++k)
                        {
                            const Vector3<T>    slavePos  = pCand + (qCand >> localPos[k]);
                            const Quaternion<T> slaveQuat = qCand * localQuat[k];
                            position[k]                   = slavePos;
                            orientation[k]                = slaveQuat;
                            LC.addParticleToCell(k, cells->computeCellHash(slavePos));
                        }
                    }
                }
            }

            GAssert(
                placed,
                std::string("Failed to place a particle without overlap after too many attempts.")
                    + " Only inserted " + std::to_string(i) + " out of "
                    + std::to_string(numParticles) + " particles. ");
            reportProgress(i);
        }
    }
    Gout("Inserted", std::to_string(numParticles), "particles.");
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class Insertion<float>;
template class Insertion<double>;