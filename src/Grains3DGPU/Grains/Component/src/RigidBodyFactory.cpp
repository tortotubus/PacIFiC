#include "RigidBodyFactory.hh"
#include "Box.hh"
#include "Cone.hh"
#include "Convex.hh"
#include "ConvexFactory.hh"
#include "Cylinder.hh"
#include "Dodecahedron.hh"
#include "Drum.hh"
#include "GrainsParameters.hh"
#include "Matrix3.hh"
#include "Octahedron.hh"
#include "QuaternionMath.hh"
#include "Rectangle.hh"
#include "RigidBody.hh"
#include "Sphere.hh"
#include "Superquadric.hh"
#include "Trapezoid.hh"
#include "Triangle.hh"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <unordered_set>
#include <vector>

/* ============================================================================================== */
/* Low-Level Methods                                                                              */
/* ============================================================================================== */
// GPU kernel to construct rigid bodies on device in batches
// Takes array of indices where each rigid body should be placed
template <typename T, typename... Arguments>
__GLOBAL__ void createRigidBodiesBatchKernel(RigidBody<T>** rb,
                                             const uint*    indices,
                                             T              crustThickness,
                                             T              density,
                                             uint           material,
                                             ConvexType     convexType,
                                             uint           count,
                                             Arguments... args)
{
    uint idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx >= count)
        return;

    Convex<T>* convex = nullptr;

    if constexpr(sizeof...(args) == 1)
    {
        if(convexType == SPHERE)
            convex = new Sphere<T>(args...);
        else if(convexType == OCTAHEDRON)
            convex = new Octahedron<T>(args...);
        else if(convexType == DODECAHEDRON)
            convex = new Dodecahedron<T>(args...);
    }
    else if constexpr(sizeof...(args) == 2)
    {
        if(convexType == CYLINDER)
            convex = new Cylinder<T>(args...);
        else if(convexType == DRUM)
            convex = new Drum<T>(args...);
        else if(convexType == CONE)
            convex = new Cone<T>(args...);
        else if(convexType == RECTANGLE)
            convex = new Rectangle<T>(args...);
        else if(convexType == TRIANGLE)
            convex = new Triangle<T>(args...);
    }
    else if constexpr(sizeof...(args) == 3)
    {
        if(convexType == BOX)
            convex = new Box<T>(args...);
        else if(convexType == TRAPEZOID)
            convex = new Trapezoid<T>(args...);
    }
    else if constexpr(sizeof...(args) == 5)
    {
        if(convexType == SUPERQUADRIC)
            convex = new Superquadric<T>(args...);
    }

    if(convex)
    {
        uint targetIdx = indices[idx];
        rb[targetIdx]  = new RigidBody<T>(convex, crustThickness, density, material);
    }
}

// -------------------------------------------------------------------------------------------------
/** @brief Fills d_RB[indices[i]] = d_RB[canonicalIdx] for all i in [0, count).
    Used to propagate a single shared device RigidBody pointer to all alias slots
    of the same type, avoiding one device allocation per particle instance. */
template <typename T>
__GLOBAL__ void fillSharedPointersKernel(RigidBody<T>** d_RB,
                                         const uint*    indices,
                                         uint           canonicalIdx,
                                         uint           count)
{
    uint i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i >= count)
        return;
    d_RB[indices[i]] = d_RB[canonicalIdx];
}

/** @brief Calls delete on every RigidBody pointer. */
template <typename T>
__GLOBAL__ void deleteRigidBodiesKernel(RigidBody<T>** d_RB, uint n)
{
    uint i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i < n)
    {
        delete d_RB[i];
        d_RB[i] = nullptr;
    }
}

/* ============================================================================================== */
/* High-Level Methods                                                                             */
/* ============================================================================================== */
// Creates and stores a RigidBody object in the host memory.
template <typename T>
__HOST__ void RigidBodyFactory<T>::create(DOMNode*         obstacles,
                                          DOMNode*         particles,
                                          ParticleData<T>& obstacleData,
                                          ParticleData<T>& particleData,
                                          uint&            numObstacles,
                                          uint&            numParticles)
{
    // Obstacles
    numObstacles              = 0;
    DOMNodeList* allObstacles = obstacles ? ReaderXML::getNodes(obstacles) : nullptr;
    uint         numRawObs    = allObstacles ? allObstacles->getLength() : 0;

    // Helper: given an <Obstacle> XML node, detect if its <Convex> child is a <ConicalWall> or
    // <CylinderWall> XML type (auto-expanded to N panels), and return N (default 32).
    // Returns 0 for standard (non-expanded) obstacles.
    auto getPanelCount = [](DOMNode* nObs, std::string& shapeType) -> uint {
        DOMNode* shape = ReaderXML::getNode(nObs, "Convex");
        if(!shape)
            return 0;
        DOMNode* shapeChild = ReaderXML::getNodeNext(shape);
        if(!shapeChild)
            return 0;
        shapeType = ReaderXML::getNodeName(shapeChild);
        if(shapeType == "ConicalWall" || shapeType == "CylinderWall")
        {
            if(ReaderXML::hasNodeAttr(shapeChild, "N"))
                return (uint)ReaderXML::getNodeAttr_Int(shapeChild, "N");
            return 32u;
        }
        return 0;
    };

    // Helper: lookup or register material, return material ID
    auto getMaterialId = [](const std::string& matStr) -> uint {
        if(GrainsParameters<T>::m_materialMap.count(matStr) == 0)
        {
            uint id = GrainsParameters<T>::m_materialMap.size();
            GrainsParameters<T>::m_materialMap.emplace(matStr, id);
        }
        return GrainsParameters<T>::m_materialMap[matStr];
    };

    // First pass: count element nodes, total sub-body slots, and build per-template flat offsets.
    // templateFlatOffset[obsIdx] = first flat component index for that template's panels.
    uint              numRefObs  = 0;
    uint              totalSlots = 0;
    std::vector<uint> templateFlatOffset;
    for(uint i = 0; i < numRawObs; ++i)
    {
        DOMNode* nObs = allObstacles->item(i);
        if(nObs->getNodeType() != 1)
            continue;
        templateFlatOffset.push_back(totalSlots);
        ++numRefObs;
        std::string shapeType;
        uint        N = getPanelCount(nObs, shapeType);
        totalSlots += (N > 0) ? N : 1u;
    }

    obstacleData.numTemplates = numRefObs;
    if(numRefObs > 0)
    {
        obstacleData.refRB.initialize(totalSlots);
        obstacleData.subBodiesPerTemplate.initialize(numRefObs);
        obstacleData.numEachRef.initialize(numRefObs);
        obstacleData.isComposite.initialize(numRefObs);
        obstacleData.refLocalPos.initialize(totalSlots);
        obstacleData.refLocalQuat.initialize(totalSlots);
        obstacleData.refInitialPos.initialize(numRefObs);
        obstacleData.refInitialOri.initialize(numRefObs);
        // Count total valid interval entries for buffer pre-allocation.
        // Expanded walls (N panels) emit N events per loading interval — one per panel.
        uint totalEvents = 0;
        for(uint i = 0; i < numRawObs; ++i)
        {
            DOMNode* nObs2 = allObstacles->item(i);
            if(nObs2->getNodeType() != 1)
                continue;
            std::string shapeType2;
            uint        panelN = getPanelCount(nObs2, shapeType2);
            uint        mult   = (panelN > 0) ? panelN : 1u;
            DOMNode*    nOL2   = ReaderXML::getNode(nObs2, "ObstacleLoading");
            if(nOL2)
            {
                DOMNodeList* ivs  = ReaderXML::getNodes(nOL2);
                uint         nRaw = ivs ? ivs->getLength() : 0;
                for(uint k = 0; k < nRaw; ++k)
                    if(ivs->item(k)->getNodeType() == 1)
                        totalEvents += mult;
            }
        }
        obstacleData.refObstacleEvents.initialize(totalEvents);

        uint obsIdx  = 0;  // per-template index
        uint slotIdx = 0;  // per-sub-body (refRB) index
        uint evIdx   = 0;
        for(uint i = 0; i < numRawObs; ++i)
        {
            DOMNode* nObs = allObstacles->item(i);
            if(nObs->getNodeType() != 1)
                continue;

            // Read world-frame transformation
            Vector3<T>    centre(T(0), T(0), T(0));
            Quaternion<T> rotation(T(0), T(0), T(0), T(1));
            DOMNode*      nTransform = ReaderXML::getNode(nObs, "Transformation");
            if(nTransform)
            {
                DOMNode* nCentre = ReaderXML::getNode(nTransform, "Centre");
                if(nCentre)
                    centre = Vector3<T>(nCentre);
                DOMNode* nRotation = ReaderXML::getNode(nTransform, "AngularPosition");
                if(nRotation)
                    rotation = Quaternion<T>(nRotation);
            }
            obstacleData.refInitialPos[obsIdx] = centre;
            obstacleData.refInitialOri[obsIdx] = rotation;
            obstacleData.numEachRef[obsIdx]    = 1;

            // Check if this is an expanded wall obstacle
            std::string shapeType;
            uint        N = getPanelCount(nObs, shapeType);

            if(N > 0)
            {
                // Expanded obstacle: create N panel sub-bodies
                DOMNode*    shape      = ReaderXML::getNode(nObs, "Convex");
                DOMNode*    shapeChild = ReaderXML::getNodeNext(shape);
                T           ct         = T(ReaderXML::getNodeAttr_Double(shape, "CrustThickness"));
                std::string matStr     = ReaderXML::getNodeAttr_String(nObs, "Material");
                uint        matId      = getMaterialId(matStr);
                const T     PI2_N      = T(2) * PI<T> / T(N);

                obstacleData.subBodiesPerTemplate[obsIdx] = N;
                obstacleData.isComposite[obsIdx]          = 1u;

                if(shapeType == "ConicalWall")
                {
                    T rB     = T(ReaderXML::getNodeAttr_Double(shapeChild, "RadiusBottom"));
                    T rT     = T(ReaderXML::getNodeAttr_Double(shapeChild, "RadiusTop"));
                    T height = T(ReaderXML::getNodeAttr_Double(shapeChild, "Height"));
                    // Panel half-widths: arc-chord = r * sin(π/N)
                    T sinPN          = sin(PI<T> / T(N));
                    T cosPN          = cos(PI<T> / T(N));
                    T wb             = rB * sinPN;  // half-width at bottom
                    T wt             = rT * sinPN;  // half-width at top
                    T bottomInradius = rB * cosPN;
                    T topInradius    = rT * cosPN;
                    T deltaInradius  = topInradius - bottomInradius;
                    T panelHeight    = sqrt(height * height + deltaInradius * deltaInradius);
                    T invPanelHeight = T(1) / panelHeight;
                    // Panel center radial distance is the average of the bottom/top chord centres.
                    T rMid = (bottomInradius + topInradius) * T(0.5);

                    for(uint k = 0; k < N; ++k)
                    {
                        T theta = PI2_N * T(k);
                        T cosT  = cos(theta);
                        T sinT  = sin(theta);

                        // Panel center in obstacle body frame: (rMid*sin(θ), 0, rMid*cos(θ))
                        // (Y is the cone axis; panel faces inward)
                        Vector3<T> localCenter(rMid * sinT, T(0), rMid * cosT);

                        // Panel orientation in obstacle body frame:
                        // Local X = tangential direction.
                        // Local Y = generator direction from the bottom chord centre to the top
                        // one. Local Z = outward normal = X x Y. Matrix3 is row-major [row0; row1;
                        // row2], columns = basis vectors.
                        Matrix3<T>    rotMat(cosT,
                                          deltaInradius * invPanelHeight * sinT,
                                          height * invPanelHeight * sinT,
                                          T(0),
                                          height * invPanelHeight,
                                          -deltaInradius * invPanelHeight,
                                          -sinT,
                                          deltaInradius * invPanelHeight * cosT,
                                          height * invPanelHeight * cosT);
                        Quaternion<T> localQuat(rotMat);

                        // Trapezoid takes full widths and the slant height between bottom and top
                        // rims.
                        Convex<T>* convex = new Trapezoid<T>(T(2) * wb, T(2) * wt, panelHeight);
                        obstacleData.refRB[slotIdx] = new RigidBody<T>(convex, ct, T(0), matId);
                        obstacleData.refLocalPos[slotIdx]  = localCenter;
                        obstacleData.refLocalQuat[slotIdx] = localQuat;
                        ++slotIdx;
                    }
                }
                else  // CylinderWall
                {
                    T radius   = T(ReaderXML::getNodeAttr_Double(shapeChild, "Radius"));
                    T height   = T(ReaderXML::getNodeAttr_Double(shapeChild, "Height"));
                    T sinPN    = sin(PI<T> / T(N));
                    T cosPN    = cos(PI<T> / T(N));
                    T panelW   = T(2) * radius * sinPN;  // full panel width
                    T inradius = radius * cosPN;         // panel center distance from axis

                    for(uint k = 0; k < N; ++k)
                    {
                        T theta = PI2_N * T(k);
                        T cosT  = cos(theta);
                        T sinT  = sin(theta);

                        Vector3<T>    localCenter(inradius * sinT, T(0), inradius * cosT);
                        Matrix3<T>    rotMat(cosT, T(0), sinT, T(0), T(1), T(0), -sinT, T(0), cosT);
                        Quaternion<T> localQuat(rotMat);

                        // Rectangle takes full dimensions: panelW, height
                        Convex<T>* convex           = new Rectangle<T>(panelW, height);
                        obstacleData.refRB[slotIdx] = new RigidBody<T>(convex, ct, T(0), matId);
                        obstacleData.refLocalPos[slotIdx]  = localCenter;
                        obstacleData.refLocalQuat[slotIdx] = localQuat;
                        ++slotIdx;
                    }
                }
            }
            else
            {
                // Standard single-body obstacle
                obstacleData.subBodiesPerTemplate[obsIdx] = 1;
                obstacleData.isComposite[obsIdx]          = 0u;
                obstacleData.refRB[slotIdx]               = new RigidBody<T>(nObs);
                obstacleData.refLocalPos[slotIdx]         = Vector3<T>(T(0), T(0), T(0));
                obstacleData.refLocalQuat[slotIdx]        = Quaternion<T>(T(0), T(0), T(0), T(1));
                ++slotIdx;
            }

            // Parse optional prescribed motion.
            // For expanded walls (N panels), each loading interval produces N events — one per
            // panel — so that all panels receive the same velocity command simultaneously.
            DOMNode* nObsLoading = ReaderXML::getNode(nObs, "ObstacleLoading");
            if(nObsLoading)
            {
                uint         firstPanelId = templateFlatOffset[obsIdx];
                uint         panelCount   = obstacleData.subBodiesPerTemplate[obsIdx];
                DOMNodeList* allIntervals = ReaderXML::getNodes(nObsLoading);
                uint         nRawIv       = allIntervals ? allIntervals->getLength() : 0;
                for(uint k = 0; k < nRawIv; ++k)
                {
                    DOMNode* nIv = allIntervals->item(k);
                    if(nIv->getNodeType() != 1)
                        continue;
                    T          tStart = T(ReaderXML::getNodeAttr_Double(nIv, "tStart"));
                    T          tEnd   = T(ReaderXML::getNodeAttr_Double(nIv, "tEnd"));
                    Vector3<T> linearVelocity(T(0), T(0), T(0));
                    Vector3<T> angularVelocity(T(0), T(0), T(0));
                    DOMNode*   nLinVel = ReaderXML::getNode(nIv, "LinearVelocity");
                    if(nLinVel)
                        linearVelocity = Vector3<T>(nLinVel);
                    DOMNode* nAngVel = ReaderXML::getNode(nIv, "AngularVelocity");
                    if(nAngVel)
                        angularVelocity = Vector3<T>(nAngVel);
                    // One event per panel
                    for(uint p = 0; p < panelCount; ++p)
                    {
                        ObstacleMotionEvent<T> ev;
                        ev.obstacleId                           = firstPanelId + p;
                        ev.tStart                               = tStart;
                        ev.tEnd                                 = tEnd;
                        ev.linearVelocity                       = linearVelocity;
                        ev.angularVelocity                      = angularVelocity;
                        obstacleData.refObstacleEvents[evIdx++] = ev;
                    }
                }
            }

            numObstacles += obstacleData.subBodiesPerTemplate[obsIdx];
            ++obsIdx;
        }
    }

    // Particles
    numParticles              = 0;
    DOMNodeList* allParticles = particles ? ReaderXML::getNodes(particles) : nullptr;
    uint         numRefPar    = allParticles ? allParticles->getLength() : 0;
    particleData.numTemplates = 0;
    if(numRefPar > 0)
    {
        // Intermediate storage (vectors before copying into GrainsMemBuffers)
        std::vector<RigidBody<T>*> rbVec;         // per-sub-body prototype
        std::vector<Vector3<T>>    localPosVec;   // per-sub-body prototype
        std::vector<Quaternion<T>> localQuatVec;  // per-sub-body prototype
        std::vector<uint>          subBodiesVec;  // per-template
        std::vector<uint>          numEachVec;    // per-template
        std::vector<uint>          isCompVec;     // per-template
        std::vector<Vector3<T>>    initPosVec;    // per-template
        std::vector<Quaternion<T>> initOriVec;    // per-template
        uint                       numTemplates = 0;

        for(uint i = 0; i < numRefPar; ++i)
        {
            DOMNode* nParticle = allParticles->item(i);
            if(nParticle->getNodeType() != 1)
                continue;
            std::string type = ReaderXML::hasNodeAttr(nParticle, "Type")
                                   ? ReaderXML::getNodeAttr_String(nParticle, "Type")
                                   : "";
            uint        M_i  = ReaderXML::getNodeAttr_Int(nParticle, "Number");

            // World transform for this template
            Vector3<T>    centre(T(0), T(0), T(0));
            Quaternion<T> rotation(T(0), T(0), T(0), T(1));
            DOMNode*      nTransform = ReaderXML::getNode(nParticle, "Transformation");
            if(nTransform)
            {
                DOMNode* nCentre = ReaderXML::getNode(nTransform, "Centre");
                if(nCentre)
                    centre = Vector3<T>(nCentre);
                DOMNode* nRotation = ReaderXML::getNode(nTransform, "AngularPosition");
                if(nRotation)
                    rotation = Quaternion<T>(nRotation);
            }

            if(type != "Composite")
            {
                // Standalone particle
                rbVec.push_back(new RigidBody<T>(nParticle));
                localPosVec.push_back(Vector3<T>(T(0), T(0), T(0)));
                localQuatVec.push_back(Quaternion<T>(T(0), T(0), T(0), T(1)));
                subBodiesVec.push_back(1u);
                numEachVec.push_back(M_i);
                isCompVec.push_back(0u);
                initPosVec.push_back(centre);
                initOriVec.push_back(rotation);
                numParticles += M_i;
                ++numTemplates;
            }
            else
            {
                // Composite particle -- collect sub-bodies
                DOMNodeList*          allChildren = ReaderXML::getNodes(nParticle);
                uint                  numChildren = allChildren ? allChildren->getLength() : 0;
                std::vector<DOMNode*> subBodies;
                for(uint k = 0; k < numChildren; ++k)
                {
                    DOMNode*    child    = allChildren->item(k);
                    std::string nodeName = ReaderXML::getNodeName(child);
                    if(nodeName == "SubBody")
                        subBodies.push_back(child);
                }
                GAssert(!subBodies.empty(), "Composite particle has no SubBody children!");

                uint S_i = subBodies.size();
                subBodiesVec.push_back(S_i);
                numEachVec.push_back(M_i);
                isCompVec.push_back(1u);
                initPosVec.push_back(centre);
                initOriVec.push_back(rotation);

                // Density and Material are read from the parent <Particle> node;
                // individual <SubBody> nodes may override Density but share Material.
                T           subDensity = ReaderXML::hasNodeAttr(nParticle, "Density")
                                             ? T(ReaderXML::getNodeAttr_Double(nParticle, "Density"))
                                             : T(0);
                std::string parentMat  = ReaderXML::getNodeAttr_String(nParticle, "Material");
                if(GrainsParameters<T>::m_materialMap.count(parentMat) == 0)
                    GrainsParameters<T>::m_materialMap.emplace(
                        parentMat,
                        (uint)GrainsParameters<T>::m_materialMap.size());
                uint matId = GrainsParameters<T>::m_materialMap[parentMat];

                for(DOMNode* nSub : subBodies)
                {
                    DOMNode*   nConvex = ReaderXML::getNode(nSub, "Convex");
                    T          ct     = T(ReaderXML::getNodeAttr_Double(nConvex, "CrustThickness"));
                    Convex<T>* convex = ConvexFactory<T>::create(nConvex);
                    T          d      = ReaderXML::hasNodeAttr(nSub, "Density")
                                            ? T(ReaderXML::getNodeAttr_Double(nSub, "Density"))
                                            : subDensity;
                    rbVec.push_back(new RigidBody<T>(convex, ct, d, matId));

                    Vector3<T>    localPos(T(0), T(0), T(0));
                    Quaternion<T> localQuat(T(0), T(0), T(0), T(1));
                    DOMNode*      nLocalTf = ReaderXML::getNode(nSub, "LocalTransformation");
                    if(nLocalTf)
                    {
                        DOMNode* nLC = ReaderXML::getNode(nLocalTf, "Centre");
                        if(nLC)
                            localPos = Vector3<T>(nLC);
                        DOMNode* nLR = ReaderXML::getNode(nLocalTf, "AngularPosition");
                        if(nLR)
                            localQuat = Quaternion<T>(nLR);
                    }
                    localPosVec.push_back(localPos);
                    localQuatVec.push_back(localQuat);
                }

                // Compute composite mass + full 3x3 inertia in master frame, then
                // diagonalize via Jacobi so the master uses principal-axis inertia.
                // All local transforms are updated to reflect the new frame; the
                // initial orientation is adjusted by the inverse rotation so world
                // placement is unchanged.
                {
                    uint base = static_cast<uint>(rbVec.size()) - S_i;

                    // Accumulate mass and 3x3 inertia tensor in master body frame
                    T totalMass = T(0);
                    T I[3][3]   = {};
                    for(uint k = 0; k < S_i; ++k)
                    {
                        RigidBody<T>* rb = rbVec[base + k];
                        T             mk = rb->getMass();
                        totalMass += mk;

                        // Rotate sub-body's diagonal inertia into master frame via q_local
                        T Ik[3];
                        rb->getInertia(Ik);
                        Matrix3<T> Rk = localQuatVec[base + k].toMatrix();
                        // I_contribution[i][j] = sum_l R(i,l) * Ik[l] * R(j,l)
                        for(int ii = 0; ii < 3; ++ii)
                            for(int jj = 0; jj < 3; ++jj)
                                for(int l = 0; l < 3; ++l)
                                    I[ii][jj] += Rk(ii, l) * Ik[l] * Rk(jj, l);

                        // Parallel-axis theorem: I += mk*(|r|^2 * delta_ij - r_i*r_j)
                        Vector3<T> rk = localPosVec[base + k];
                        T          r2 = rk[X] * rk[X] + rk[Y] * rk[Y] + rk[Z] * rk[Z];
                        I[0][0] += mk * (r2 - rk[X] * rk[X]);
                        I[1][1] += mk * (r2 - rk[Y] * rk[Y]);
                        I[2][2] += mk * (r2 - rk[Z] * rk[Z]);
                        I[0][1] -= mk * rk[X] * rk[Y];
                        I[1][0] = I[0][1];
                        I[0][2] -= mk * rk[X] * rk[Z];
                        I[2][0] = I[0][2];
                        I[1][2] -= mk * rk[Y] * rk[Z];
                        I[2][1] = I[1][2];
                    }

                    // Jacobi eigendecomposition (cyclic sweeps, max 50 sweeps).
                    // V accumulates Givens rotations; columns of V are eigenvectors.
                    // After convergence, V^T * I * V = diag(eigenvalues).
                    T V[3][3] = {{T(1), T(0), T(0)}, {T(0), T(1), T(0)}, {T(0), T(0), T(1)}};
                    for(int sweep = 0; sweep < 50; ++sweep)
                    {
                        bool converged = true;
                        for(int p = 0; p < 2; ++p)
                        {
                            for(int q = p + 1; q < 3; ++q)
                            {
                                T Ipq = I[p][q];
                                if(std::abs(Ipq) < T(1e-14))
                                    continue;
                                converged = false;
                                T tau     = (I[p][p] - I[q][q]) / (T(2) * Ipq);
                                T t = (tau >= T(0)) ? T(1) / (tau + std::sqrt(T(1) + tau * tau))
                                                    : T(1) / (tau - std::sqrt(T(1) + tau * tau));
                                T c = T(1) / std::sqrt(T(1) + t * t);
                                T s = t * c;
                                // Update diagonal (Ipq -> 0)
                                I[p][p] -= t * Ipq;
                                I[q][q] += t * Ipq;
                                I[p][q] = I[q][p] = T(0);
                                // Update remaining off-diagonal row/col
                                for(int r = 0; r < 3; ++r)
                                {
                                    if(r == p || r == q)
                                        continue;
                                    T Irp   = c * I[r][p] - s * I[r][q];
                                    T Irq   = s * I[r][p] + c * I[r][q];
                                    I[r][p] = I[p][r] = Irp;
                                    I[r][q] = I[q][r] = Irq;
                                }
                                // Accumulate eigenvector matrix: V_new = V_old * G
                                for(int r = 0; r < 3; ++r)
                                {
                                    T Vrp   = c * V[r][p] - s * V[r][q];
                                    T Vrq   = s * V[r][p] + c * V[r][q];
                                    V[r][p] = Vrp;
                                    V[r][q] = Vrq;
                                }
                            }
                        }
                        if(converged)
                            break;
                    }

                    // q_p rotates from old master frame to principal frame (= V^T).
                    // Matrix3 layout is row-major: buf[3*i+j] = M(i,j).
                    T             Rp_buf[9] = {V[0][0],
                                               V[1][0],
                                               V[2][0],
                                               V[0][1],
                                               V[1][1],
                                               V[2][1],
                                               V[0][2],
                                               V[1][2],
                                               V[2][2]};
                    Matrix3<T>    Rp(Rp_buf);
                    Quaternion<T> q_p(Rp);
                    T             qn = norm(q_p);
                    if(qn > T(1e-12))
                        q_p *= (T(1) / qn);

                    // Rotate local positions and local quaternions into principal frame
                    for(uint k = 0; k < S_i; ++k)
                    {
                        localPosVec[base + k]  = q_p >> localPosVec[base + k];
                        localQuatVec[base + k] = q_p * localQuatVec[base + k];
                        T qnk                  = norm(localQuatVec[base + k]);
                        if(qnk > T(1e-12))
                            localQuatVec[base + k] *= (T(1) / qnk);
                    }

                    // Compensate initial orientation: q_init_new = q_init_old * conj(q_p)
                    // so the composite's world placement is unchanged.
                    const Quaternion<T>& cq_p = q_p;
                    initOriVec.back()         = initOriVec.back() * conjugate(cq_p);
                    T qni                     = norm(initOriVec.back());
                    if(qni > T(1e-12))
                        initOriVec.back() *= (T(1) / qni);

                    // Store composite mass + principal inertia in the master prototype (k=0).
                    // All M_i instances are copied from this prototype, so they all
                    // get the correct properties without repeating the computation.
                    rbVec[base]->setCompositeProperties(totalMass, I[0][0], I[1][1], I[2][2]);
                }

                numParticles += M_i * S_i;
                ++numTemplates;
            }
        }

        // Copy from vectors into GrainsMemBuffers
        uint totalSubBodyProtos   = rbVec.size();
        particleData.numTemplates = numTemplates;
        if(numTemplates > 0)
        {
            particleData.subBodiesPerTemplate.initialize(numTemplates);
            particleData.numEachRef.initialize(numTemplates);
            particleData.isComposite.initialize(numTemplates);
            particleData.refInitialPos.initialize(numTemplates);
            particleData.refInitialOri.initialize(numTemplates);
        }
        if(totalSubBodyProtos > 0)
        {
            particleData.refRB.initialize(totalSubBodyProtos);
            particleData.refLocalPos.initialize(totalSubBodyProtos);
            particleData.refLocalQuat.initialize(totalSubBodyProtos);
        }
        for(uint t = 0; t < numTemplates; ++t)
        {
            particleData.subBodiesPerTemplate[t] = subBodiesVec[t];
            particleData.numEachRef[t]           = numEachVec[t];
            particleData.isComposite[t]          = isCompVec[t];
            particleData.refInitialPos[t]        = initPosVec[t];
            particleData.refInitialOri[t]        = initOriVec[t];
        }
        for(uint r = 0; r < totalSubBodyProtos; ++r)
        {
            particleData.refRB[r]        = rbVec[r];
            particleData.refLocalPos[r]  = localPosVec[r];
            particleData.refLocalQuat[r] = localQuatVec[r];
        }
    }
}

// -------------------------------------------------------------------------------------------------
// Constructs RigidBody objects on device by grouping similar bodies and launching one kernel per
// group
template <typename T>
__HOST__ void
    RigidBodyFactory<T>::copyHostToDevice(GrainsMemBuffer<RigidBody<T>*, MemType::HOST>&   h_RB,
                                          GrainsMemBuffer<RigidBody<T>*, MemType::DEVICE>& d_RB)
{
    d_RB.initialize(h_RB.getSize());

    uint numRB = h_RB.getSize();
    if(numRB == 0)
        return;

    // Define a structure to identify unique rigid body types
    struct RBProperties
    {
        ConvexType type;
        T          density;
        T          crustThickness;
        uint       material;
        // Parameters for convex shapes (max 5 for superquadric)
        T    params[5];
        uint numParams;

        bool operator==(const RBProperties& other) const
        {
            if(type != other.type || material != other.material || numParams != other.numParams)
                return false;
            if(std::abs(density - other.density) > T(1e-9)
               || std::abs(crustThickness - other.crustThickness) > T(1e-9))
                return false;
            for(uint i = 0; i < numParams; ++i)
                if(std::abs(params[i] - other.params[i]) > T(1e-9))
                    return false;
            return true;
        }
    };

    // Step 1: Analyze h_RB to find unique types and their indices
    std::vector<RBProperties>      uniqueTypes;
    std::vector<std::vector<uint>> indicesPerType;

    for(uint i = 0; i < numRB; ++i)
    {
        RBProperties prop{};
        Convex<T>*   convex = h_RB[i]->getConvex();
        prop.type           = convex->getConvexType();

        auto snapshot       = h_RB[i]->getPropertiesSnapshot();
        prop.density        = snapshot.mass / snapshot.convex->computeVolume();
        prop.crustThickness = snapshot.crustThickness;
        prop.material       = snapshot.material;

        // Extract parameters based on type
        prop.numParams = 0;
        if(prop.type == SPHERE)
        {
            Sphere<T>* s   = dynamic_cast<Sphere<T>*>(convex);
            prop.params[0] = s->getRadius();
            prop.numParams = 1;
        }
        else if(prop.type == OCTAHEDRON)
        {
            Octahedron<T>* o = dynamic_cast<Octahedron<T>*>(convex);
            prop.params[0]   = o->getRadius();
            prop.numParams   = 1;
        }
        else if(prop.type == DODECAHEDRON)
        {
            Dodecahedron<T>* d = dynamic_cast<Dodecahedron<T>*>(convex);
            prop.params[0]     = d->getRadius();
            prop.numParams     = 1;
        }
        else if(prop.type == BOX)
        {
            Box<T>*    b   = dynamic_cast<Box<T>*>(convex);
            Vector3<T> L   = b->getExtent();
            prop.params[0] = L[X];
            prop.params[1] = L[Y];
            prop.params[2] = L[Z];
            prop.numParams = 3;
        }
        else if(prop.type == CYLINDER)
        {
            Cylinder<T>* c = dynamic_cast<Cylinder<T>*>(convex);
            prop.params[0] = c->getRadius();
            prop.params[1] = c->getHeight();
            prop.numParams = 2;
        }
        else if(prop.type == DRUM)
        {
            Drum<T>* d     = dynamic_cast<Drum<T>*>(convex);
            prop.params[0] = d->getRadius();
            prop.params[1] = d->getHeight();
            prop.numParams = 2;
        }
        else if(prop.type == CONE)
        {
            Cone<T>* c     = dynamic_cast<Cone<T>*>(convex);
            prop.params[0] = c->getRadius();
            prop.params[1] = c->getHeight();
            prop.numParams = 2;
        }
        else if(prop.type == RECTANGLE)
        {
            Rectangle<T>* r = dynamic_cast<Rectangle<T>*>(convex);
            Vector3<T>    L = r->getExtent();  // returns half-extents (m_LX, m_LY)
            prop.params[0]  = T(2) * L[X];     // pass full extents; Rectangle ctor halves them
            prop.params[1]  = T(2) * L[Y];
            prop.numParams  = 2;
        }
        else if(prop.type == TRAPEZOID)
        {
            Trapezoid<T>* t = dynamic_cast<Trapezoid<T>*>(convex);
            prop.params[0]  = T(2) * t->getHalfWidthBottom();
            prop.params[1]  = T(2) * t->getHalfWidthTop();
            prop.params[2]  = T(2) * t->getHalfHeight();
            prop.numParams  = 3;
        }
        else if(prop.type == TRIANGLE)
        {
            Triangle<T>* t = dynamic_cast<Triangle<T>*>(convex);
            prop.params[0] = T(2) * t->getHalfBase();
            prop.params[1] = T(2) * t->getHalfHeight();
            prop.numParams = 2;
        }
        else if(prop.type == SUPERQUADRIC)
        {
            Superquadric<T>* sq = dynamic_cast<Superquadric<T>*>(convex);
            Vector3<T>       L  = sq->getExtent();
            Vector3<T>       N  = sq->getExponent();
            prop.params[0]      = L[X];
            prop.params[1]      = L[Y];
            prop.params[2]      = L[Z];
            prop.params[3]      = N[X];
            prop.params[4]      = N[Y];
            prop.numParams      = 5;
        }

        // Find if this type already exists
        int typeIdx = -1;
        for(uint j = 0; j < uniqueTypes.size(); ++j)
        {
            if(uniqueTypes[j] == prop)
            {
                typeIdx = j;
                break;
            }
        }

        // Add to existing type or create new type
        if(typeIdx >= 0)
        {
            indicesPerType[typeIdx].push_back(i);
        }
        else
        {
            uniqueTypes.push_back(prop);
            indicesPerType.push_back({i});
        }
    }

    // Set device heap size to accommodate UNIQUE types only (not one entry per particle).
    // All particles of the same type share a single device RigidBody via pointer aliasing.
    size_t requiredHeap = uniqueTypes.size() * 1024;  // 1KB per unique type
    size_t currentHeap;
    cudaDeviceGetLimit(&currentHeap, cudaLimitMallocHeapSize);
    if(requiredHeap > currentHeap)
    {
        cudaDeviceSetLimit(cudaLimitMallocHeapSize, requiredHeap);
    }

    // Step 2: For each unique type, create exactly ONE device RigidBody (at the canonical
    // index = first occurrence), then propagate that pointer to all alias slots of the same
    // type using fillSharedPointersKernel.  This replaces the old per-particle device-new
    // loop and reduces kernel launches from O(N/2048) to O(numUniqueTypes).
    for(uint typeIdx = 0; typeIdx < uniqueTypes.size(); ++typeIdx)
    {
        const RBProperties&      prop    = uniqueTypes[typeIdx];
        const std::vector<uint>& indices = indicesPerType[typeIdx];
        uint                     count   = indices.size();

        // Allocate device memory for all indices of this type
        uint* d_indices;
        cudaErrCheck(cudaMalloc(&d_indices, count * sizeof(uint)));
        cudaErrCheck(
            cudaMemcpy(d_indices, indices.data(), count * sizeof(uint), cudaMemcpyHostToDevice));

        // Create ONE device RigidBody at the canonical slot (indices[0]) using a 1-thread kernel
        {
            uint numThreads = 1;
            uint numBlocks  = 1;

            // Launch kernel based on type
            if(prop.type == SPHERE)
            {
                createRigidBodiesBatchKernel<<<numBlocks, numThreads>>>(d_RB.getData(),
                                                                        d_indices,
                                                                        prop.crustThickness,
                                                                        prop.density,
                                                                        prop.material,
                                                                        SPHERE,
                                                                        1u,
                                                                        prop.params[0]);
            }
            else if(prop.type == OCTAHEDRON)
            {
                createRigidBodiesBatchKernel<<<numBlocks, numThreads>>>(d_RB.getData(),
                                                                        d_indices,
                                                                        prop.crustThickness,
                                                                        prop.density,
                                                                        prop.material,
                                                                        OCTAHEDRON,
                                                                        1u,
                                                                        prop.params[0]);
            }
            else if(prop.type == DODECAHEDRON)
            {
                createRigidBodiesBatchKernel<<<numBlocks, numThreads>>>(d_RB.getData(),
                                                                        d_indices,
                                                                        prop.crustThickness,
                                                                        prop.density,
                                                                        prop.material,
                                                                        DODECAHEDRON,
                                                                        1u,
                                                                        prop.params[0]);
            }
            else if(prop.type == BOX)
            {
                createRigidBodiesBatchKernel<<<numBlocks, numThreads>>>(d_RB.getData(),
                                                                        d_indices,
                                                                        prop.crustThickness,
                                                                        prop.density,
                                                                        prop.material,
                                                                        BOX,
                                                                        1u,
                                                                        prop.params[0],
                                                                        prop.params[1],
                                                                        prop.params[2]);
            }
            else if(prop.type == CYLINDER)
            {
                createRigidBodiesBatchKernel<<<numBlocks, numThreads>>>(d_RB.getData(),
                                                                        d_indices,
                                                                        prop.crustThickness,
                                                                        prop.density,
                                                                        prop.material,
                                                                        CYLINDER,
                                                                        1u,
                                                                        prop.params[0],
                                                                        prop.params[1]);
            }
            else if(prop.type == DRUM)
            {
                createRigidBodiesBatchKernel<<<numBlocks, numThreads>>>(d_RB.getData(),
                                                                        d_indices,
                                                                        prop.crustThickness,
                                                                        prop.density,
                                                                        prop.material,
                                                                        DRUM,
                                                                        1u,
                                                                        prop.params[0],
                                                                        prop.params[1]);
            }
            else if(prop.type == CONE)
            {
                createRigidBodiesBatchKernel<<<numBlocks, numThreads>>>(d_RB.getData(),
                                                                        d_indices,
                                                                        prop.crustThickness,
                                                                        prop.density,
                                                                        prop.material,
                                                                        CONE,
                                                                        1u,
                                                                        prop.params[0],
                                                                        prop.params[1]);
            }
            else if(prop.type == RECTANGLE)
            {
                createRigidBodiesBatchKernel<<<numBlocks, numThreads>>>(d_RB.getData(),
                                                                        d_indices,
                                                                        prop.crustThickness,
                                                                        prop.density,
                                                                        prop.material,
                                                                        RECTANGLE,
                                                                        1u,
                                                                        prop.params[0],
                                                                        prop.params[1]);
            }
            else if(prop.type == TRAPEZOID)
            {
                createRigidBodiesBatchKernel<<<numBlocks, numThreads>>>(d_RB.getData(),
                                                                        d_indices,
                                                                        prop.crustThickness,
                                                                        prop.density,
                                                                        prop.material,
                                                                        TRAPEZOID,
                                                                        1u,
                                                                        prop.params[0],
                                                                        prop.params[1],
                                                                        prop.params[2]);
            }
            else if(prop.type == TRIANGLE)
            {
                createRigidBodiesBatchKernel<<<numBlocks, numThreads>>>(d_RB.getData(),
                                                                        d_indices,
                                                                        prop.crustThickness,
                                                                        prop.density,
                                                                        prop.material,
                                                                        TRIANGLE,
                                                                        1u,
                                                                        prop.params[0],
                                                                        prop.params[1]);
            }
            else if(prop.type == SUPERQUADRIC)
            {
                createRigidBodiesBatchKernel<<<numBlocks, numThreads>>>(d_RB.getData(),
                                                                        d_indices,
                                                                        prop.crustThickness,
                                                                        prop.density,
                                                                        prop.material,
                                                                        SUPERQUADRIC,
                                                                        1u,
                                                                        prop.params[0],
                                                                        prop.params[1],
                                                                        prop.params[2],
                                                                        prop.params[3],
                                                                        prop.params[4]);
            }
            else
            {
                cudaFree(d_indices);
                GAbort("Convex type is not implemented for GPU! Aborting Grains!");
            }

        }

        cudaErrCheck(cudaDeviceSynchronize());
        cudaError_t err = cudaGetLastError();
        if(err != cudaSuccess)
        {
            cudaFree(d_indices);
            GAbort(("CUDA error in RigidBody creation: "
                    + std::string(cudaGetErrorString(err)))
                       .c_str());
        }

        // Propagate the canonical pointer (d_RB[indices[0]]) to all alias slots
        uint aliases = count - 1;
        if(aliases > 0)
        {
            uint canonicalIdx = indices[0];
            uint numThreads   = 256;
            uint numBlocks    = (aliases + numThreads - 1) / numThreads;
            fillSharedPointersKernel<<<numBlocks, numThreads>>>(
                d_RB.getData(), d_indices + 1, canonicalIdx, aliases);
            cudaErrCheck(cudaDeviceSynchronize());
        }

        cudaFree(d_indices);
    }

    cudaErrCheck(cudaDeviceSynchronize());
}

// -------------------------------------------------------------------------------------------------
template <typename T>
__HOST__ void RigidBodyFactory<T>::freeDevice(GrainsMemBuffer<RigidBody<T>*, MemType::DEVICE>& d_RB)
{
    const uint n = static_cast<uint>(d_RB.getSize());
    if(n == 0)
        return;

    // copyHostToDevice stores shared pointers: many slots in d_RB point to the SAME device
    // object (one per unique particle type).  Copy all pointer values to the host, deduplicate,
    // then delete only the canonical (unique) objects to avoid double-free.
    std::vector<RigidBody<T>*> h_ptrs(n);
    cudaErrCheck(cudaMemcpy(
        h_ptrs.data(), d_RB.getData(), n * sizeof(RigidBody<T>*), cudaMemcpyDeviceToHost));

    // Deduplicate: collect unique non-null pointers
    std::unordered_set<RigidBody<T>*> seen;
    std::vector<RigidBody<T>*>        unique_ptrs;
    for(auto* p : h_ptrs)
        if(p != nullptr && seen.insert(p).second)
            unique_ptrs.push_back(p);

    if(unique_ptrs.empty())
        return;

    // Upload the deduplicated list to device and run the delete kernel
    const uint numUnique = static_cast<uint>(unique_ptrs.size());
    RigidBody<T>** d_unique;
    cudaErrCheck(cudaMalloc(&d_unique, numUnique * sizeof(RigidBody<T>*)));
    cudaErrCheck(cudaMemcpy(
        d_unique, unique_ptrs.data(), numUnique * sizeof(RigidBody<T>*), cudaMemcpyHostToDevice));

    const uint numThreads = 256;
    const uint numBlocks  = (numUnique + numThreads - 1) / numThreads;
    deleteRigidBodiesKernel<<<numBlocks, numThreads>>>(d_unique, numUnique);
    cudaErrCheck(cudaDeviceSynchronize());
    cudaFree(d_unique);
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class RigidBodyFactory<float>;
template class RigidBodyFactory<double>;
