#include "ComponentManager.hh"
#include "BodyTag.hh"
#include "ComponentManagerCommon.hh"
#include "ComponentManagerGPU_Kernels.hh"
#include "ForceModule.hh"
#include "QuaternionMath.hh"
#include "VectorMath.hh"

// -------------------------------------------------------------------------------------------------
template <typename T, MemType M>
ComponentManager<T, M>::ComponentManager(
    GrainsMemBuffer<RigidBody<T>*, M>*                            rigidBody,
    GrainsMemBuffer<uint, MemType::HOST>&&                        bodyTags,
    GrainsMemBuffer<Vector3<T>, MemType::HOST>&&                  position,
    GrainsMemBuffer<Quaternion<T>, MemType::HOST>&&               orientation,
    GrainsMemBuffer<Vector3<T>, MemType::HOST>&&                  localPos,
    GrainsMemBuffer<Quaternion<T>, MemType::HOST>&&               localQuat,
    const GrainsMemBuffer<ObstacleMotionEvent<T>, MemType::HOST>& obstacleEvents,
    uint                                                          nObstacles,
    uint                                                          nParticles,
    uint                                                          nComposites,
    uint                                                          nSubBodies)
    : m_rigidBody(rigidBody)
    , m_velocity(nParticles + nObstacles)
    , m_torce(nParticles + nObstacles)
    , m_masterSlot(nComposites + 1)
    , m_compositeVelocity(nComposites + 1)
    , m_compositeTorce(nComposites + 1)
    , m_counts{nObstacles, nParticles, 0u, nComposites, nSubBodies}
{
    const uint nComp = nParticles + nObstacles;
    GAssert(rigidBody->getSize() == nComp, "Rigid body size mismatch");
    const uint* hbtPtr = bodyTags.getData();  // capture before potential move

    GrainsParameters<T>::m_hasMovingObstacles = (obstacleEvents.getSize() > 0);

    if constexpr(M == MemType::HOST)
    {
        // pointer transfer
        m_bodyTag    = std::move(bodyTags);
        m_position   = std::move(position);
        m_quaternion = std::move(orientation);
        m_localPos   = std::move(localPos);
        m_localQuat  = std::move(localQuat);

        // Populate masterSlot (m_bodyTag is HOST-accessible)
        for(uint i = 0; i < nComp; ++i)
        {
            const uint tag = m_bodyTag[i];
            if(isSubBody(tag) && getSubBodyLocalIdx(tag) == 0u)
                m_masterSlot[getCompositeIdx(tag)] = i;
        }

        // Initialise composite frame state from master slot data
        m_compositePosition.initialize(nComposites + 1);
        m_compositeQuaternion.initialize(nComposites + 1);
        initCompositeFrames();

        // Derive sub-body world transforms and velocities from composite frame
        updateSubBodyPositions();
    }
    else
    {
        // Upload all HOST buffers to device
        m_bodyTag.initialize(nComp);
        m_bodyTag.copyFrom(bodyTags);
        m_position.initialize(nComp);
        m_position.copyFrom(position);
        m_quaternion.initialize(nComp);
        m_quaternion.copyFrom(orientation);
        m_localPos.initialize(nComp);
        m_localPos.copyFrom(localPos);
        m_localQuat.initialize(nComp);
        m_localQuat.copyFrom(localQuat);

        // Build masterSlot on host (bodyTags still HOST-accessible; not moved), then upload
        GrainsMemBuffer<uint, MemType::HOST>          hMasterSlot(nComposites + 1);
        GrainsMemBuffer<Vector3<T>, MemType::HOST>    hCompositePosition(nComposites + 1);
        GrainsMemBuffer<Quaternion<T>, MemType::HOST> hCompositeQuaternion(nComposites + 1);
        for(uint i = 0; i < nComp; ++i)
        {
            const uint tag = bodyTags[i];
            if(isSubBody(tag) && getSubBodyLocalIdx(tag) == 0u)
            {
                const uint cIdx            = getCompositeIdx(tag);
                hMasterSlot[cIdx]          = i;
                hCompositePosition[cIdx]   = position[i];     // position[mSlot] == CM at insertion
                hCompositeQuaternion[cIdx] = orientation[i];  // matching CM quaternion
            }
        }
        m_masterSlot.copyFrom(hMasterSlot);
        m_compositePosition.initialize(nComposites + 1);
        m_compositePosition.copyFrom(hCompositePosition);
        m_compositeQuaternion.initialize(nComposites + 1);
        m_compositeQuaternion.copyFrom(hCompositeQuaternion);

        cudaErrCheck(cudaStreamCreateWithFlags(&m_obstacleMoveStream, cudaStreamNonBlocking));
    }

    // Upload the flat event buffer
    if(obstacleEvents.getSize() > 0)
    {
        m_obstacleEvents.initialize(obstacleEvents.getSize());
        m_obstacleEvents.copyFrom(obstacleEvents);
        if constexpr(M == MemType::DEVICE)
        {
            m_obstaclesMovedFlag.initialize(1);
            *m_obstaclesMovedFlag.getData() = 0;
        }
    }

    // Create CollisionDetectionModule
    m_collisionDetectionModule = std::make_unique<CollisionDetectionModule<T, M>>(
        m_rigidBody,
        m_position,
        m_quaternion,
        m_bodyTag,
        hbtPtr,
        GrainsParameters<T>::m_collisionDetection,
        m_counts.numObstacles,
        m_counts.numParticles);

    // Size pair-indexed buffers to CDM's initial pair capacity
    const size_t pairCapacity = m_collisionDetectionModule->getPairBufferSize();
    m_contactInfo.initialize(pairCapacity);
    m_pairList.initialize(pairCapacity);

    // Create ForceModule
    m_forceModule = std::make_unique<ForceModule<T, M>>(pairCapacity,
                                                        GrainsParameters<T>::m_isContactWithMemory);
}

// -------------------------------------------------------------------------------------------------
template <typename T, MemType M>
ComponentManager<T, M>::~ComponentManager()
{
    if constexpr(M == MemType::DEVICE)
    {
        if(m_obstacleMoveStream != nullptr)
            cudaErrCheck(cudaStreamDestroy(m_obstacleMoveStream));
    }
}

// -------------------------------------------------------------------------------------------------
template <typename T, MemType M>
const GrainsMemBuffer<Vector3<T>, M>& ComponentManager<T, M>::getLocalPos() const
{
    return m_localPos;
}

// -------------------------------------------------------------------------------------------------
template <typename T, MemType M>
const GrainsMemBuffer<Quaternion<T>, M>& ComponentManager<T, M>::getLocalQuat() const
{
    return m_localQuat;
}

// -------------------------------------------------------------------------------------------------
template <typename T, MemType M>
const GrainsMemBuffer<uint, M>& ComponentManager<T, M>::getBodyTag() const
{
    return m_bodyTag;
}

// -------------------------------------------------------------------------------------------------
template <typename T, MemType M>
const GrainsMemBuffer<Vector3<T>, M>& ComponentManager<T, M>::getPosition() const
{
    return m_position;
}

// -------------------------------------------------------------------------------------------------
template <typename T, MemType M>
const GrainsMemBuffer<Quaternion<T>, M>& ComponentManager<T, M>::getQuaternion() const
{
    return m_quaternion;
}

// -------------------------------------------------------------------------------------------------
template <typename T, MemType M>
const GrainsMemBuffer<Kinematics<T>, M>& ComponentManager<T, M>::getVelocity() const
{
    return m_velocity;
}

// -------------------------------------------------------------------------------------------------
template <typename T, MemType M>
const GrainsMemBuffer<Torce<T>, M>& ComponentManager<T, M>::getTorce() const
{
    return m_torce;
}

// -------------------------------------------------------------------------------------------------
template <typename T, MemType M>
const NeighborList<T, M>* ComponentManager<T, M>::getNeighborList() const
{
    return m_collisionDetectionModule ? m_collisionDetectionModule->getNeighborList() : nullptr;
}

// -------------------------------------------------------------------------------------------------
template <typename T, MemType M>
const CollisionDetectionModule<T, M>* ComponentManager<T, M>::getCollisionDetectionModule() const
{
    return m_collisionDetectionModule.get();
}

// -------------------------------------------------------------------------------------------------
template <typename T, MemType M>
uint ComponentManager<T, M>::getNumberOfParticles() const
{
    return m_counts.numParticles;
}

// -------------------------------------------------------------------------------------------------
template <typename T, MemType M>
uint ComponentManager<T, M>::getNumberOfObstacles() const
{
    return m_counts.numObstacles;
}

// -------------------------------------------------------------------------------------------------
template <typename T, MemType M>
uint ComponentManager<T, M>::getNumberOfComposites() const
{
    return m_counts.numComposites;
}

// -------------------------------------------------------------------------------------------------
template <typename T, MemType M>
uint ComponentManager<T, M>::getNumberOfSubBodies() const
{
    return m_counts.numSubBodies;
}

// -------------------------------------------------------------------------------------------------
template <typename T, MemType M>
uint ComponentManager<T, M>::getNumberOfPairs() const
{
    return m_counts.numPairs;
}

// -------------------------------------------------------------------------------------------------
template <typename T, MemType M>
void ComponentManager<T, M>::copyTo_PostProcessing(
    const std::unique_ptr<ComponentManager<T, MemType::HOST>>& other)
{
    other->setPosition(m_position);
    other->setQuaternion(m_quaternion);
    other->setVelocity(m_velocity);
}

// -------------------------------------------------------------------------------------------------
// Snapshots compositePosition, compositeQuaternion, and compositeVelocity from the master
// slot's current position, quaternion, and velocity for every composite.
// Call once after construction/insertion (HOST only).
template <typename T, MemType M>
void ComponentManager<T, M>::initCompositeFrames()
{
    if constexpr(M == MemType::HOST)
    {
        if(m_counts.numComposites == 0)
            return;
        const uint nComp = m_counts.numObstacles + m_counts.numParticles;
        for(uint i = 0; i < nComp; ++i)
        {
            const uint tag = m_bodyTag[i];
            if(isSubBody(tag) && getSubBodyLocalIdx(tag) == 0u)
            {
                const uint cIdx             = getCompositeIdx(tag);
                m_compositePosition[cIdx]   = m_position[i];
                m_compositeQuaternion[cIdx] = m_quaternion[i];
                m_compositeVelocity[cIdx]   = m_velocity[i];
            }
        }
    }
}

// -------------------------------------------------------------------------------------------------
template <typename T, MemType M>
void ComponentManager<T, M>::insertParticles(const std::unique_ptr<Insertion<T>>& insertionPolicy)
{
    // Insertion is a HOST-only operation. The DEVICE branch is discarded by if constexpr so
    // that type-incompatible calls to insertionPolicy->insert() are never compiled for DEVICE.
    if constexpr(M == MemType::HOST)
    {
        insertionPolicy->insert(m_rigidBody,
                                m_position,
                                m_quaternion,
                                m_velocity,
                                GrainsParameters<T>::m_collisionDetection.linkedCellParameters,
                                m_counts.numObstacles,
                                m_counts.numParticles,
                                m_bodyTag,
                                m_localPos,
                                m_localQuat);
        // Re-read composite frame origins: after array / force-insertion, position[masterSlot]
        // holds the inserted composite-center (CM) positions; we snapshot them here so the
        // subsequent updateSubBodyPositions() call can derive every sub-body world transform
        // correctly.
        initCompositeFrames();
    }
}

// -------------------------------------------------------------------------------------------------
template <typename T, MemType M>
void ComponentManager<T, M>::detectCollisions()
{
    m_collisionDetectionModule->run(*m_rigidBody,
                                    m_position,
                                    m_quaternion,
                                    m_velocity,
                                    m_torce,
                                    m_bodyTag,
                                    m_localPos,
                                    m_localQuat,
                                    m_masterSlot,
                                    m_contactInfo,
                                    m_pairList,
                                    m_counts);
    // Keep ForceModule GPU compaction buffers (m_activeIndex, m_cubSelectTempStorage)
    // in sync with the pair buffer that the CDModule may have just grown.
    m_forceModule->resizeBuffers(m_pairList.getSize());
}

// -------------------------------------------------------------------------------------------------
template <typename T, MemType M>
void ComponentManager<T, M>::computeContactForces(
    const GrainsMemBuffer<ContactForceModel<T>*, M>& CF)
{
    m_forceModule->run(CF,
                       m_rigidBody,
                       m_position,
                       m_compositePosition,
                       m_velocity,
                       m_pairList,
                       m_contactInfo,
                       m_torce,
                       m_compositeTorce,
                       m_bodyTag,
                       m_counts);
}
// --------------------------------------------------------------------------------------------------
// Advances all moving components (particles + obstacles) by one time step.
template <typename T, MemType M>
void ComponentManager<T, M>::moveComponents(const GrainsMemBuffer<TimeIntegrator<T>*, M>& TI)
{
    using GP = GrainsParameters<T>;

    const T cellSize = GP::m_collisionDetection.linkedCellParameters.minCellSize
                       * GP::m_collisionDetection.linkedCellParameters.cellSizeFactor;
    const T    maxDisplacementSq = cellSize * cellSize;
    const T    dt                = GP::m_dt;
    const T    currentTime       = GP::m_simulationState.time - dt;
    const bool noRelinkObstacles = GP::m_collisionDetection.linkedCellParameters.noRelinkObstacles;
    GP::m_simulationState.obstaclesMoved = false;

    if constexpr(M == MemType::HOST)
    {
        // Obstacles
        if(GP::m_hasMovingObstacles)
        {
            for(uint obstacleId = 0; obstacleId < m_counts.numObstacles; ++obstacleId)
                m_velocity[obstacleId] = Kinematics<T>();

            for(uint i = 0; i < m_obstacleEvents.getSize(); ++i)
            {
                if(moveObstacle_common(m_position.getData(),
                                       m_quaternion.getData(),
                                       m_velocity.getData(),
                                       m_obstacleEvents[i],
                                       currentTime,
                                       dt))
                    GP::m_simulationState.obstaclesMoved = true;
            }
        }

        // Particles
        for(uint pID = m_counts.numObstacles; pID < m_counts.numObstacles + m_counts.numParticles;
            ++pID)
        {
            if(isSubBody(m_bodyTag[pID]) && getSubBodyLocalIdx(m_bodyTag[pID]) != 0u)
                continue;
            moveParticles_common(TI.getData(),
                                 m_rigidBody->getData(),
                                 m_position.getData(),
                                 m_compositePosition.getData(),
                                 m_quaternion.getData(),
                                 m_compositeQuaternion.getData(),
                                 m_velocity.getData(),
                                 m_compositeVelocity.getData(),
                                 m_torce.getData(),
                                 m_compositeTorce.getData(),
                                 m_bodyTag.getData(),
                                 pID,
                                 maxDisplacementSq);
        }
        updateSubBodyPositions();
    }
    else  // DEVICE
    {
        // Obstacles
        if(GP::m_hasMovingObstacles)
        {
            const uint numEvents            = static_cast<uint>(m_obstacleEvents.getSize());
            *m_obstaclesMovedFlag.getData() = 0;

            if(m_counts.numObstacles > 0)
            {
                cudaErrCheck(cudaMemsetAsync(m_velocity.getData(),
                                             0,
                                             m_counts.numObstacles * sizeof(Kinematics<T>),
                                             m_obstacleMoveStream));
            }

            // Launch thread-per-event kernel
            moveObstacles_Kernel<<<1, numEvents, 0, m_obstacleMoveStream>>>(
                m_obstacleEvents.getData(),
                numEvents,
                currentTime,
                dt,
                m_position.getData(),
                m_quaternion.getData(),
                m_velocity.getData(),
                m_obstaclesMovedFlag.getDeviceData());
        }

        // Particles
        uint numThreads, numBlocks;
        computeOptimalThreadsAndBlocks(m_counts.numParticles, GP::m_GPU, numBlocks, numThreads);
        moveParticles_Kernel<<<numBlocks, numThreads>>>(TI.getData(),
                                                        m_rigidBody->getData(),
                                                        m_position.getData(),
                                                        m_compositePosition.getData(),
                                                        m_quaternion.getData(),
                                                        m_compositeQuaternion.getData(),
                                                        m_velocity.getData(),
                                                        m_compositeVelocity.getData(),
                                                        m_torce.getData(),
                                                        m_compositeTorce.getData(),
                                                        m_bodyTag.getData(),
                                                        m_counts.numObstacles,
                                                        m_counts.numParticles,
                                                        maxDisplacementSq);
        updateSubBodyPositions();

        cudaErrCheck(cudaStreamSynchronize(0));
        if(GP::m_hasMovingObstacles)
        {
            cudaErrCheck(cudaStreamSynchronize(m_obstacleMoveStream));
            GP::m_simulationState.obstaclesMoved = (*m_obstaclesMovedFlag.getData() != 0);
        }
    }

    if(noRelinkObstacles)
        GP::m_simulationState.obstaclesMoved = false;
}

// -------------------------------------------------------------------------------------------------
// Performs the second velocity half-kick (KDK Step 3; no-op for single-pass schemes).
template <typename T, MemType M>
void ComponentManager<T, M>::advanceVelocity(const GrainsMemBuffer<TimeIntegrator<T>*, M>& TI)
{
    if constexpr(M == MemType::HOST)
    {
        for(uint pID = m_counts.numObstacles; pID < m_counts.numObstacles + m_counts.numParticles;
            ++pID)
        {
            if(isSubBody(m_bodyTag[pID]) && getSubBodyLocalIdx(m_bodyTag[pID]) != 0u)
                continue;
            advanceVelocity_common(TI.getData(),
                                   m_rigidBody->getData(),
                                   m_quaternion.getData(),
                                   m_compositeQuaternion.getData(),
                                   m_velocity.getData(),
                                   m_compositeVelocity.getData(),
                                   m_torce.getData(),
                                   m_compositeTorce.getData(),
                                   m_bodyTag.getData(),
                                   pID);
        }
    }
    else
    {
        uint numThreads, numBlocks;
        computeOptimalThreadsAndBlocks(m_counts.numParticles,
                                       GrainsParameters<T>::m_GPU,
                                       numBlocks,
                                       numThreads);
        advanceVelocity_Kernel<<<numBlocks, numThreads>>>(TI.getData(),
                                                          m_rigidBody->getData(),
                                                          m_quaternion.getData(),
                                                          m_compositeQuaternion.getData(),
                                                          m_velocity.getData(),
                                                          m_compositeVelocity.getData(),
                                                          m_torce.getData(),
                                                          m_compositeTorce.getData(),
                                                          m_bodyTag.getData(),
                                                          m_counts.numObstacles,
                                                          m_counts.numParticles);
    }
    // Re-derive sub-body velocities so m_velocity[k] reflects the updated compositeVelocity
    updateSubBodyPositions();
}

// -------------------------------------------------------------------------------------------------
// Slaves non-master sub-body world transforms to their composite master.
template <typename T, MemType M>
void ComponentManager<T, M>::updateSubBodyPositions()
{
    if(m_counts.numSubBodies == 0)
        return;

    if constexpr(M == MemType::HOST)
    {
        const uint nTotal = m_counts.numObstacles + m_counts.numParticles;
        for(uint cID = m_counts.numObstacles; cID < nTotal; ++cID)
        {
            const uint tag = m_bodyTag[cID];
            if(!isSubBody(tag))
                continue;
            const uint           cIdx = getCompositeIdx(tag);
            const Vector3<T>&    cp   = m_compositePosition[cIdx];
            const Quaternion<T>& cq   = m_compositeQuaternion[cIdx];
            m_position[cID]           = cp + (cq >> m_localPos[cID]);
            m_quaternion[cID]         = cq * m_localQuat[cID];
            T qn                      = norm(m_quaternion[cID]);
            if(qn > T(1e-12))
                m_quaternion[cID] *= (T(1) / qn);
            // Derive per-sub-body velocity: v_k = v_CM + omega x r_k
            const Vector3<T>&    r     = m_position[cID] - cp;
            const Kinematics<T>& compK = m_compositeVelocity[cIdx];
            m_velocity[cID]            = Kinematics<T>(compK.getTranslationalComponent()
                                                + (compK.getAngularComponent() ^ r),
                                            compK.getAngularComponent());
        }
    }
    else
    {
        const uint nTotal = m_counts.numObstacles + m_counts.numParticles;
        uint       numThreads, numBlocks;
        computeOptimalThreadsAndBlocks(nTotal, GrainsParameters<T>::m_GPU, numBlocks, numThreads);
        updateSubBodies_Kernel<<<numBlocks, numThreads>>>(m_localPos.getData(),
                                                          m_localQuat.getData(),
                                                          m_position.getData(),
                                                          m_compositePosition.getData(),
                                                          m_quaternion.getData(),
                                                          m_compositeQuaternion.getData(),
                                                          m_velocity.getData(),
                                                          m_compositeVelocity.getData(),
                                                          m_bodyTag.getData(),
                                                          nTotal);
    }
}

// --------------------------------------------------------------------------------------------------
// Explicit instantiations
template class ComponentManager<float, MemType::HOST>;
template class ComponentManager<double, MemType::HOST>;
template class ComponentManager<float, MemType::DEVICE>;
template class ComponentManager<double, MemType::DEVICE>;
