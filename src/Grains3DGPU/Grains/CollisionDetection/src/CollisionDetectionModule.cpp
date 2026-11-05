#include "GrainsPlatform.hh"
#if GRAINS_CUDA_COMPILER
#include <cub/cub.cuh>
#endif
#include <vector>

#include "BodyTag.hh"
#include "CollisionDetectionCommon.hh"
#include "CollisionComputation.hh"
#include "CollisionComputationCPU.hh"
#include "CollisionComputationCUDA.hh"
#include "CollisionDetectionModule.hh"
#if GRAINS_CUDA_COMPILER
#include "CollisionDetectionModule_Kernels.hh"
#endif
#include "GrainsParameters.hh"
#include "GrainsUtils.hh"

// -------------------------------------------------------------------------------------------------
// Constructor: builds the NeighborList object, then allocates pair-indexed buffers
template <typename T, MemType M>
CollisionDetectionModule<T, M>::CollisionDetectionModule(
    const GrainsMemBuffer<RigidBody<T>*, M>* rigidBody,
    const GrainsMemBuffer<Vector3<T>, M>&    positions,
    const GrainsMemBuffer<Quaternion<T>, M>& orientations,
    const GrainsMemBuffer<uint, M>&          bodyTags,
    const uint*                              hostBodyTags,
    const CollisionDetectionParameters<T>&   CD,
    uint                                     nObstacles,
    uint                                     nParticles)
    : m_particleSorter(nObstacles, nParticles)
    , m_neighborList(NeighborListFactory<T, M>::create(
          rigidBody, positions, orientations, CD, nObstacles, nParticles))
{
    size_t freeMem;
    if constexpr(M == MemType::HOST)
        freeMem = getAvailableHostMemory();
#if GRAINS_CUDA_COMPILER
    else
        freeMem = getAvailableDeviceMemory();
#endif

    const auto& LCD       = GrainsParameters<T>::m_collisionDetection.linkedCellParameters;
    size_t estimatedPairs = static_cast<size_t>(nParticles) * LCD.initialNumberOfPairsPerParticle;
    size_t maxPairs       = static_cast<size_t>(nObstacles) * nParticles
                      + static_cast<size_t>(nParticles) * (nParticles - 1) / 2;
    estimatedPairs = std::min(estimatedPairs, maxPairs);

    size_t sizePerPair = sizeof(Vector3<T>) + sizeof(Quaternion<T>) + sizeof(ContactInfo<T>)
                         + sizeof(uint8_t) + sizeof(uint);
    size_t sizeNeeded = estimatedPairs * sizePerPair;
    GAssert(sizeNeeded < freeMem,
            "Not enough memory for pair-dependent buffers in CollisionDetectionModule!");

    m_relPosition.initialize(estimatedPairs);
    m_relQuaternion.initialize(estimatedPairs);
    m_contactInfoLocal.initialize(estimatedPairs);
#if GRAINS_CUDA_COMPILER
    if constexpr(M == MemType::DEVICE)
    {
        m_bvPassFlags.initialize(estimatedPairs);
        m_bvPassPairIndices.initialize(estimatedPairs);
        m_bvPassPairCountMapped.initialize(1);
        *m_bvPassPairCountMapped.getData() = 0;
        // Query CUB scratch size for DeviceSelect::Flagged over estimatedPairs elements
        cub::CountingInputIterator<uint> countIter(0);
        cub::DeviceSelect::Flagged(nullptr,
                                   m_cubTempStorageBytes,
                                   countIter,
                                   (uint8_t*)nullptr,
                                   (uint*)nullptr,
                                   (int*)nullptr,
                                   (int)estimatedPairs);
        m_cubTempStorage.initialize(m_cubTempStorageBytes);
    }
#endif

    // Always build the compact shape/BV tables so that toggling usePrebuiltShapes at
    // runtime (e.g. across benchmark iterations) never reads an uninitialized m_shapeData.
    // Cost is O(nUniqueShapes) -- dominated by every other init step.
    {
        const uint nComponents = nObstacles + nParticles;
        buildShapeAndBVData(rigidBody->getData(), hostBodyTags, nComponents);
    }
}

// -------------------------------------------------------------------------------------------------
// Gets the current pair list pointer from the neighbor list
template <typename T, MemType M>
const IndexPair* CollisionDetectionModule<T, M>::getPairList() const
{
    return m_neighborList->getData();
}

// -------------------------------------------------------------------------------------------------
// Gets the current number of active pairs in the neighbor list
template <typename T, MemType M>
uint CollisionDetectionModule<T, M>::getPairCount() const
{
    return m_neighborList->getSize();
}

// -------------------------------------------------------------------------------------------------
// Gets the allocated size of the pair buffers (may exceed getPairCount())
template <typename T, MemType M>
size_t CollisionDetectionModule<T, M>::getPairBufferSize() const
{
    return m_relPosition.getSize();
}

// -------------------------------------------------------------------------------------------------
// Gets the neighbor list (read-only)
template <typename T, MemType M>
const NeighborList<T, M>* CollisionDetectionModule<T, M>::getNeighborList() const
{
    return m_neighborList.get();
}

// -------------------------------------------------------------------------------------------------
// Runs the full collision detection pipeline
template <typename T, MemType M>
void CollisionDetectionModule<T, M>::run(const RigidBody<T>* const*          rigidBodies,
                                         GrainsMemBuffer<Vector3<T>, M>&     positions,
                                         GrainsMemBuffer<Quaternion<T>, M>&  orientations,
                                         GrainsMemBuffer<Kinematics<T>, M>&  velocities,
                                         GrainsMemBuffer<Torce<T>, M>&       torces,
                                         GrainsMemBuffer<uint, M>&           bodyTags,
                                         GrainsMemBuffer<Vector3<T>, M>&     localPos,
                                         GrainsMemBuffer<Quaternion<T>, M>&  localQuat,
                                         GrainsMemBuffer<uint, M>&           masterSlot,
                                         GrainsMemBuffer<ContactInfo<T>, M>& contactInfo,
                                         GrainsMemBuffer<IndexPair, M>&      pairList,
                                         ComponentCounts&                    counts)
{
    auto& gt = GrainsParameters<T>::m_cdmTimer;
    gt.start(CDMStage::Total);
    sortParticles(positions,
                  orientations,
                  velocities,
                  torces,
                  bodyTags,
                  localPos,
                  localQuat,
                  masterSlot,
                  counts);
    updateNeighborList(positions, pairList, contactInfo, counts);
    if(GrainsParameters<T>::m_collisionDetection.useRelativeTransformations)
    {
        computeRelativeTransformations(positions, orientations, pairList);
        detectCollisionsComponents(rigidBodies, pairList, bodyTags, counts);
        transformContactInfo(positions, orientations, pairList, contactInfo);
    }
    else
        detectCollisionsComponentsGlobal(rigidBodies,
                                         positions,
                                         orientations,
                                         pairList,
                                         contactInfo,
                                         bodyTags,
                                         counts);
    gt.stop(CDMStage::Total);
}

// -------------------------------------------------------------------------------------------------
// Resizes per-pair buffers and ComponentManager's contactInfo / pairList
template <typename T, MemType M>
void CollisionDetectionModule<T, M>::resizePairBuffers(
    GrainsMemBuffer<IndexPair, M>&      pairList,
    GrainsMemBuffer<ContactInfo<T>, M>& contactInfo,
    uint                                size)
{
    m_relPosition.resize(size);
    m_relQuaternion.resize(size);
    m_contactInfoLocal.resize(size);
#if GRAINS_CUDA_COMPILER
    if constexpr(M == MemType::DEVICE)
    {
        m_bvPassFlags.resize(size);
        m_bvPassPairIndices.resize(size);
        // Re-query CUB scratch size for the new capacity
        cub::CountingInputIterator<uint> countIter(0);
        size_t                           newBytes = 0;
        cub::DeviceSelect::Flagged(nullptr,
                                   newBytes,
                                   countIter,
                                   (uint8_t*)nullptr,
                                   (uint*)nullptr,
                                   (int*)nullptr,
                                   (int)size);
        if(newBytes > m_cubTempStorageBytes)
        {
            m_cubTempStorageBytes = newBytes;
            m_cubTempStorage.resize(m_cubTempStorageBytes);
        }
    }
#endif

    contactInfo.resize(size);
    pairList.resize(size);
}

// -------------------------------------------------------------------------------------------------
// Updates the neighbor list; resizes buffers and increments global counter when rebuilt
template <typename T, MemType M>
void CollisionDetectionModule<T, M>::updateNeighborList(
    GrainsMemBuffer<Vector3<T>, M>&     positions,
    GrainsMemBuffer<IndexPair, M>&      pairList,
    GrainsMemBuffer<ContactInfo<T>, M>& contactInfo,
    ComponentCounts&                    counts)
{
    const uint nObstacles = counts.numObstacles;
    const uint nParticles = counts.numParticles;
    uint&      numPairs   = counts.numPairs;
    auto&      SS         = GrainsParameters<T>::m_simulationState;
    auto&      gt         = GrainsParameters<T>::m_cdmTimer;

    gt.start(CDMStage::NeighborList);

    bool updated = m_neighborList->updateNeighborList(positions, nObstacles, nParticles);
    if(updated)
    {
        const uint newSize = m_neighborList->getSize();
        resizePairBuffers(pairList, contactInfo, newSize);
        m_neighborList->getBuffer().copyTo(pairList);
        numPairs = newSize;
        SS.neighborListUpdateCount++;
    }

    gt.stop(CDMStage::NeighborList);
}

// -------------------------------------------------------------------------------------------------
// Computes per-pair relative position / quaternion (B in A-local frame)
template <typename T, MemType M>
void CollisionDetectionModule<T, M>::computeRelativeTransformations(
    const GrainsMemBuffer<Vector3<T>, M>&    positions,
    const GrainsMemBuffer<Quaternion<T>, M>& orientations,
    const GrainsMemBuffer<IndexPair, M>&     pairList)
{
    auto& timer = GrainsParameters<T>::m_cdmTimer;
    timer.start(CDMStage::RelativeTransform);
    const RelativeTransformView<T> view{pairList.getData(),
                                        positions.getData(),
                                        orientations.getData(),
                                        m_relPosition.getData(),
                                        m_relQuaternion.getData(),
                                        m_neighborList->getSize()};
    if constexpr(M == MemType::HOST)
        grains::cpu::computeRelativeTransforms(view);
    else
        grains::cuda::computeRelativeTransforms(view);
    timer.stop(CDMStage::RelativeTransform);
}

// -------------------------------------------------------------------------------------------------
// Builds compact ShapeData and BVData tables indexed by shapeId (not component slot).
// Scans bodyTagsData to find nUniqueShapes and the representative slot for each shapeId.
// DEVICE: uploads repSlots to GPU, launches two fill kernels (one thread per unique shape).
// HOST:   calls buildShapeAndBVData (CPU loop over unique shapes).
template <typename T, MemType M>
void CollisionDetectionModule<T, M>::buildShapeAndBVData(const RigidBody<T>* const* rigidBodies,
                                                         const uint*                bodyTagsData,
                                                         uint                       nComponents)
{
    // Pass 1: scan bodyTags on CPU (bodyTagsData is always host-accessible here)
    constexpr uint   MAX_SHAPES = 1024u;  // 10-bit shapeId field
    std::vector<int> repSlotsV(MAX_SHAPES, -1);
    uint             nUniqueShapes = 0;
    for(uint i = 0; i < nComponents; ++i)
    {
        const uint shapeId = getShapeId(bodyTagsData[i]);
        if(repSlotsV[shapeId] < 0)
        {
            repSlotsV[shapeId] = static_cast<int>(i);
            if(shapeId + 1 > nUniqueShapes)
                nUniqueShapes = shapeId + 1;
        }
    }
    m_nUniqueShapes = nUniqueShapes;

    // Build compact repSlots array (size = nUniqueShapes)
    std::vector<uint> repSlots(nUniqueShapes);
    for(uint k = 0; k < nUniqueShapes; ++k)
        repSlots[k] = (repSlotsV[k] >= 0) ? static_cast<uint>(repSlotsV[k]) : 0u;

    // Pass 2: fill tables
    m_shapeData.initialize(nUniqueShapes);
    m_bvData.initialize(nUniqueShapes);

#if GRAINS_CUDA_COMPILER
    if constexpr(M == MemType::DEVICE)
    {
        // Upload repSlots to device
        GrainsMemBuffer<uint, MemType::DEVICE> d_repSlots;
        d_repSlots.initialize(nUniqueShapes);
        cudaMemcpy(d_repSlots.getData(),
                   repSlots.data(),
                   nUniqueShapes * sizeof(uint),
                   cudaMemcpyHostToDevice);

        uint numThreads, numBlocks;
        computeOptimalThreadsAndBlocks(nUniqueShapes,
                                       GrainsParameters<T>::m_GPU,
                                       numBlocks,
                                       numThreads);
        fillShapeData_Kernel<T><<<numBlocks, numThreads>>>(m_shapeData.getData(),
                                                           rigidBodies,
                                                           d_repSlots.getData(),
                                                           nUniqueShapes);
        fillBVData_Kernel<T><<<numBlocks, numThreads>>>(m_bvData.getData(),
                                                        rigidBodies,
                                                        d_repSlots.getData(),
                                                        nUniqueShapes);
        cudaDeviceSynchronize();
    }
    else
#endif

    {
        ::buildShapeAndBVData(m_shapeData.getData(),
                              m_bvData.getData(),
                              rigidBodies,
                              repSlots.data(),
                              nUniqueShapes);
    }
}

// -------------------------------------------------------------------------------------------------
// Performs a narrow-phase GJK-based collision detection for each pair and writes contact info in
// A-local frame. DEVICE BV-ON and composite paths use the compacted m_bvPassPairIndices list.
template <typename T, MemType M>
void CollisionDetectionModule<T, M>::detectCollisionsComponents(
    const RigidBody<T>* const*           rigidBodies,
    const GrainsMemBuffer<IndexPair, M>& pairList,
    const GrainsMemBuffer<uint, M>&      bodyTags,
    const ComponentCounts&               counts)
{
    using GP           = GrainsParameters<T>;
    const auto& config = GP::m_collisionDetection;
    auto&       timer  = GP::m_cdmTimer;
    uint        count  = m_neighborList->getSize();
    const uint* active = nullptr;
#if GRAINS_CUDA_COMPILER
    if constexpr(M == MemType::DEVICE)
    {
        if(config.boundingVolumeType != BoundingVolumeType::OFF || counts.numComposites > 0)
        {
            filterPairsBV(rigidBodies,
                          pairList,
                          m_contactInfoLocal,
                          bodyTags.getData(),
                          counts.numComposites);
            count  = static_cast<uint>(*m_bvPassPairCountMapped.getData());
            active = m_bvPassPairIndices.getData();
        }
    }
#endif
    const CollisionPairView<T>     pairs{pairList.getData(),
                                     rigidBodies,
                                     m_shapeData.getData(),
                                     bodyTags.getData(),
                                     m_contactInfoLocal.getData(),
                                     active,
                                     count,
                                     counts.numComposites};
    const RelativeCollisionView<T> view{pairs, m_relPosition.getData(), m_relQuaternion.getData()};
    NarrowPhaseOptions             options{config.narrowPhaseType,
                               config.boundingVolumeType,
                               config.gjkAcceleration,
                               config.usePrebuiltShapes};
    timer.start(CDMStage::NarrowPhase);
    if constexpr(M == MemType::HOST)
        grains::cpu::detectRelativeCollisions(view, options);
    else
    {
        // CUDA bounding-volume rejection has already run in filterPairsBV.
        options.boundingVolume = BoundingVolumeType::OFF;
        grains::cuda::detectRelativeCollisions(view, options);
    }
    timer.stop(CDMStage::NarrowPhase);
}

// -------------------------------------------------------------------------------------------------
// Narrow-phase GJK using world-frame positions/quaternions directly.
template <typename T, MemType M>
void CollisionDetectionModule<T, M>::detectCollisionsComponentsGlobal(
    const RigidBody<T>* const*               rigidBodies,
    const GrainsMemBuffer<Vector3<T>, M>&    positions,
    const GrainsMemBuffer<Quaternion<T>, M>& orientations,
    const GrainsMemBuffer<IndexPair, M>&     pairList,
    GrainsMemBuffer<ContactInfo<T>, M>&      contactInfo,
    const GrainsMemBuffer<uint, M>&          bodyTags,
    const ComponentCounts&                   counts)
{
    using GP           = GrainsParameters<T>;
    const auto& config = GP::m_collisionDetection;
    auto&       timer  = GP::m_cdmTimer;
    uint        count  = m_neighborList->getSize();
    const uint* active = nullptr;
#if GRAINS_CUDA_COMPILER
    if constexpr(M == MemType::DEVICE)
    {
        if(config.boundingVolumeType != BoundingVolumeType::OFF || counts.numComposites > 0)
        {
            if(config.useRelativeTransformations)
            {
                if(config.boundingVolumeType != BoundingVolumeType::OFF)
                    computeRelativeTransformations(positions, orientations, pairList);
                filterPairsBV(rigidBodies,
                              pairList,
                              contactInfo,
                              bodyTags.getData(),
                              counts.numComposites);
            }
            else
                filterPairsBV_global(rigidBodies,
                                     positions,
                                     orientations,
                                     pairList,
                                     contactInfo,
                                     bodyTags.getData(),
                                     counts.numComposites);
            count  = static_cast<uint>(*m_bvPassPairCountMapped.getData());
            active = m_bvPassPairIndices.getData();
        }
    }
#endif
    const CollisionPairView<T>   pairs{pairList.getData(),
                                     rigidBodies,
                                     m_shapeData.getData(),
                                     bodyTags.getData(),
                                     contactInfo.getData(),
                                     active,
                                     count,
                                     counts.numComposites};
    const GlobalCollisionView<T> view{pairs, positions.getData(), orientations.getData()};
    NarrowPhaseOptions           options{config.narrowPhaseType,
                               config.boundingVolumeType,
                               config.gjkAcceleration,
                               config.usePrebuiltShapes};
    timer.start(CDMStage::NarrowPhase);
    if constexpr(M == MemType::HOST)
        grains::cpu::detectGlobalCollisions(view, options);
    else
    {
        options.boundingVolume = BoundingVolumeType::OFF;
        grains::cuda::detectGlobalCollisions(view, options);
    }
    timer.stop(CDMStage::NarrowPhase);
}

// -------------------------------------------------------------------------------------------------
// Launches a BV kernel to write flags, then uses CUB DeviceSelect::Flagged to compact passing
// pair indices into m_bvPassPairIndices. Also serves as the composite-only filter when BV is
// OFF. Only invoked on the DEVICE path.
template <typename T, MemType M>
void CollisionDetectionModule<T, M>::filterPairsBV(
    const RigidBody<T>* const*           rigidBodies,
    const GrainsMemBuffer<IndexPair, M>& pairList,
    GrainsMemBuffer<ContactInfo<T>, M>&  contactInfoLocal,
    const uint*                          bodyTags,
    uint                                 numComposites)
{
#if GRAINS_CUDA_COMPILER
    if constexpr(M == MemType::DEVICE)
    {
        auto& gt = GrainsParameters<T>::m_cdmTimer;
        gt.start(CDMStage::BVFilter);

        const uint nPairs = m_neighborList->getSize();
        if(nPairs == 0)
        {
            *m_bvPassPairCountMapped.getData() = 0;
            gt.stop(CDMStage::BVFilter);
            return;
        }
        uint numThreads, numBlocks;
        computeOptimalThreadsAndBlocks(nPairs, GrainsParameters<T>::m_GPU, numBlocks, numThreads);

        // Step 1: BV pass/fail flags + no-contact sentinels for rejected pairs
        const BoundingVolumeType bvType
            = GrainsParameters<T>::m_collisionDetection.boundingVolumeType;
        const bool usePrebuilt = GrainsParameters<T>::m_collisionDetection.usePrebuiltShapes;
        if(bvType == BoundingVolumeType::OBC)
        {
            if(usePrebuilt)
                filterPairsBV_Kernel<T, BoundingVolumeType::OBC>
                    <<<numBlocks, numThreads>>>(m_bvData.getData(),
                                                pairList.getData(),
                                                bodyTags,
                                                numComposites,
                                                m_relPosition.getData(),
                                                m_relQuaternion.getData(),
                                                contactInfoLocal.getData(),
                                                m_bvPassFlags.getData(),
                                                nPairs);
            else
                filterPairsBV_Kernel<T, BoundingVolumeType::OBC>
                    <<<numBlocks, numThreads>>>(rigidBodies,
                                                pairList.getData(),
                                                bodyTags,
                                                numComposites,
                                                m_relPosition.getData(),
                                                m_relQuaternion.getData(),
                                                contactInfoLocal.getData(),
                                                m_bvPassFlags.getData(),
                                                nPairs);
        }
        else if(bvType == BoundingVolumeType::OBB)
        {
            if(usePrebuilt)
                filterPairsBV_Kernel<T, BoundingVolumeType::OBB>
                    <<<numBlocks, numThreads>>>(m_bvData.getData(),
                                                pairList.getData(),
                                                bodyTags,
                                                numComposites,
                                                m_relPosition.getData(),
                                                m_relQuaternion.getData(),
                                                contactInfoLocal.getData(),
                                                m_bvPassFlags.getData(),
                                                nPairs);
            else
                filterPairsBV_Kernel<T, BoundingVolumeType::OBB>
                    <<<numBlocks, numThreads>>>(rigidBodies,
                                                pairList.getData(),
                                                bodyTags,
                                                numComposites,
                                                m_relPosition.getData(),
                                                m_relQuaternion.getData(),
                                                contactInfoLocal.getData(),
                                                m_bvPassFlags.getData(),
                                                nPairs);
        }
        else
        {
            // BV is OFF: composite-only filter -- relPosition/relQuaternion not used by the
            // kernel. Prebuilt is identical for OFF, use regular kernel.
            filterPairsBV_Kernel<T, BoundingVolumeType::OFF>
                <<<numBlocks, numThreads>>>(rigidBodies,
                                            pairList.getData(),
                                            bodyTags,
                                            numComposites,
                                            nullptr,
                                            nullptr,
                                            contactInfoLocal.getData(),
                                            m_bvPassFlags.getData(),
                                            nPairs);
        }
        cudaDeviceSynchronize();

        // Step 2: CUB compaction - select original pair indices where flag == 1
        cub::CountingInputIterator<uint> countIter(0);
        cub::DeviceSelect::Flagged(m_cubTempStorage.getData(),
                                   m_cubTempStorageBytes,
                                   countIter,
                                   m_bvPassFlags.getData(),
                                   m_bvPassPairIndices.getData(),
                                   m_bvPassPairCountMapped.getDeviceData(),
                                   (int)nPairs);
        // CUB wrote the count directly into mapped pinned memory via the device alias;
        // synchronize ensures the host-side read below sees the completed value.
        cudaDeviceSynchronize();

        gt.stop(CDMStage::BVFilter);
    }
#endif
}

// -------------------------------------------------------------------------------------------------
// BV filter using world-frame positions/quaternions; computes relative transforms on-the-fly
// inside filterPairsBV_Global_Kernel to avoid a separate pre-pass.
template <typename T, MemType M>
void CollisionDetectionModule<T, M>::filterPairsBV_global(
    const RigidBody<T>* const*               rigidBodies,
    const GrainsMemBuffer<Vector3<T>, M>&    positions,
    const GrainsMemBuffer<Quaternion<T>, M>& orientations,
    const GrainsMemBuffer<IndexPair, M>&     pairList,
    GrainsMemBuffer<ContactInfo<T>, M>&      contactInfoLocal,
    const uint*                              bodyTags,
    uint                                     numComposites)
{
#if GRAINS_CUDA_COMPILER
    if constexpr(M == MemType::DEVICE)
    {
        auto& gt = GrainsParameters<T>::m_cdmTimer;
        gt.start(CDMStage::BVFilter);

        const uint nPairs = m_neighborList->getSize();
        if(nPairs == 0)
        {
            *m_bvPassPairCountMapped.getData() = 0;
            gt.stop(CDMStage::BVFilter);
            return;
        }
        uint numThreads, numBlocks;
        computeOptimalThreadsAndBlocks(nPairs, GrainsParameters<T>::m_GPU, numBlocks, numThreads);

        const BoundingVolumeType bvType
            = GrainsParameters<T>::m_collisionDetection.boundingVolumeType;
        const bool usePrebuilt = GrainsParameters<T>::m_collisionDetection.usePrebuiltShapes;
        if(bvType == BoundingVolumeType::OBC)
        {
            if(usePrebuilt)
                filterPairsBV_Global_Kernel<T, BoundingVolumeType::OBC>
                    <<<numBlocks, numThreads>>>(m_bvData.getData(),
                                                pairList.getData(),
                                                bodyTags,
                                                numComposites,
                                                positions.getData(),
                                                orientations.getData(),
                                                contactInfoLocal.getData(),
                                                m_bvPassFlags.getData(),
                                                nPairs);
            else
                filterPairsBV_Global_Kernel<T, BoundingVolumeType::OBC>
                    <<<numBlocks, numThreads>>>(rigidBodies,
                                                pairList.getData(),
                                                bodyTags,
                                                numComposites,
                                                positions.getData(),
                                                orientations.getData(),
                                                contactInfoLocal.getData(),
                                                m_bvPassFlags.getData(),
                                                nPairs);
        }
        else if(bvType == BoundingVolumeType::OBB)
        {
            if(usePrebuilt)
                filterPairsBV_Global_Kernel<T, BoundingVolumeType::OBB>
                    <<<numBlocks, numThreads>>>(m_bvData.getData(),
                                                pairList.getData(),
                                                bodyTags,
                                                numComposites,
                                                positions.getData(),
                                                orientations.getData(),
                                                contactInfoLocal.getData(),
                                                m_bvPassFlags.getData(),
                                                nPairs);
            else
                filterPairsBV_Global_Kernel<T, BoundingVolumeType::OBB>
                    <<<numBlocks, numThreads>>>(rigidBodies,
                                                pairList.getData(),
                                                bodyTags,
                                                numComposites,
                                                positions.getData(),
                                                orientations.getData(),
                                                contactInfoLocal.getData(),
                                                m_bvPassFlags.getData(),
                                                nPairs);
        }
        else
        {
            filterPairsBV_Global_Kernel<T, BoundingVolumeType::OFF>
                <<<numBlocks, numThreads>>>(rigidBodies,
                                            pairList.getData(),
                                            bodyTags,
                                            numComposites,
                                            positions.getData(),
                                            orientations.getData(),
                                            contactInfoLocal.getData(),
                                            m_bvPassFlags.getData(),
                                            nPairs);
        }
        cudaDeviceSynchronize();

        cub::CountingInputIterator<uint> countIter(0);
        cub::DeviceSelect::Flagged(m_cubTempStorage.getData(),
                                   m_cubTempStorageBytes,
                                   countIter,
                                   m_bvPassFlags.getData(),
                                   m_bvPassPairIndices.getData(),
                                   m_bvPassPairCountMapped.getDeviceData(),
                                   (int)nPairs);
        cudaDeviceSynchronize();

        gt.stop(CDMStage::BVFilter);
    }
#endif
}
template <typename T, MemType M>
void CollisionDetectionModule<T, M>::transformContactInfo(
    const GrainsMemBuffer<Vector3<T>, M>&    positions,
    const GrainsMemBuffer<Quaternion<T>, M>& orientations,
    const GrainsMemBuffer<IndexPair, M>&     pairList,
    GrainsMemBuffer<ContactInfo<T>, M>&      contactInfo)
{
    auto& timer = GrainsParameters<T>::m_cdmTimer;
    timer.start(CDMStage::Transform);
    const ContactTransformView<T> view{pairList.getData(),
                                       positions.getData(),
                                       orientations.getData(),
                                       m_contactInfoLocal.getData(),
                                       contactInfo.getData(),
                                       m_neighborList->getSize()};
    if constexpr(M == MemType::HOST)
        grains::cpu::transformContacts(view);
    else
        grains::cuda::transformContacts(view);
    timer.stop(CDMStage::Transform);
}

// -------------------------------------------------------------------------------------------------
// Sorts particles by Morton codes
template <typename T, MemType M>
void CollisionDetectionModule<T, M>::sortParticles(GrainsMemBuffer<Vector3<T>, M>&    positions,
                                                   GrainsMemBuffer<Quaternion<T>, M>& orientations,
                                                   GrainsMemBuffer<Kinematics<T>, M>& velocities,
                                                   GrainsMemBuffer<Torce<T>, M>&      torces,
                                                   GrainsMemBuffer<uint, M>&          bodyTags,
                                                   GrainsMemBuffer<Vector3<T>, M>&    localPos,
                                                   GrainsMemBuffer<Quaternion<T>, M>& localQuat,
                                                   GrainsMemBuffer<uint, M>&          masterSlot,
                                                   const ComponentCounts&             counts)
{
    const uint numComposites = counts.numComposites;
    const uint nObstacles    = counts.numObstacles;
    const uint nParticles    = counts.numParticles;
    using GP                 = GrainsParameters<T>;
    auto& SS                 = GP::m_simulationState;
    auto& LC                 = GP::m_collisionDetection.linkedCellParameters;
    auto& gt                 = GP::m_cdmTimer;

    if(LC.sortFrequency > 0 && SS.neighborListUpdateCount % LC.sortFrequency == 0)
    {
        gt.start(CDMStage::Sort);
        m_particleSorter.sortParticles(positions,
                                       velocities,
                                       orientations,
                                       torces,
                                       bodyTags,
                                       localPos,
                                       localQuat,
                                       nObstacles,
                                       nParticles);
        SS.particlesSorted = true;
        // Rebuild master-slot lookup so composite queries stay valid after the reorder
        if(numComposites > 0)
        {
            const uint nTotal = nObstacles + nParticles;
            if constexpr(M == MemType::HOST)
            {
                for(uint cID = 0; cID < nTotal; ++cID)
                {
                    const uint tag = bodyTags[cID];
                    if(isSubBody(tag) && getSubBodyLocalIdx(tag) == 0u)
                        masterSlot[getCompositeIdx(tag)] = cID;
                }
            }
#if GRAINS_CUDA_COMPILER
            else if constexpr(M == MemType::DEVICE)
            {
                uint numThreads, numBlocks;
                computeOptimalThreadsAndBlocks(nTotal,
                                               GrainsParameters<T>::m_GPU,
                                               numBlocks,
                                               numThreads);
                rebuildMasterSlot_Kernel<<<numBlocks, numThreads>>>(masterSlot.getData(),
                                                                    bodyTags.getData(),
                                                                    nTotal);
            }
#endif
        }
        gt.stop(CDMStage::Sort);
    }
    else
        SS.particlesSorted = false;
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiations
template class CollisionDetectionModule<float, MemType::HOST>;
template class CollisionDetectionModule<double, MemType::HOST>;
#if GRAINS_CUDA_COMPILER
template class CollisionDetectionModule<float, MemType::DEVICE>;
#endif
#if GRAINS_CUDA_COMPILER
template class CollisionDetectionModule<double, MemType::DEVICE>;
#endif
