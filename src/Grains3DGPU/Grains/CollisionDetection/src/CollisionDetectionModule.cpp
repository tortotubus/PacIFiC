#include <cub/cub.cuh>
#include <type_traits>
#include <vector>

#include "BodyTag.hh"
#include "CollisionDetectionCommon.hh"
#include "CollisionDetectionModule.hh"
#include "CollisionDetectionModule_Kernels.hh"
#include "GrainsParameters.hh"
#include "GrainsUtils.hh"

// -------------------------------------------------------------------------------------------------
// Dispatches from runtime (NarrowPhaseType, bool) to compile-time template arguments by calling
//   f( std::integral_constant<GJKType, GJKV>{}, std::bool_constant<GJKA>{} )
// with the four concrete instantiations. The generic lambda in the call site then recovers
// GJKV / GJKA as constexpr values via "constexpr auto V = decltype(tag)::value".
template <typename Func>
static void dispatchGJK(NarrowPhaseType npType, bool gjkAcc, Func&& f)
{
    using J  = std::integral_constant<GJKType, GJKType::JOHNSON>;
    using SV = std::integral_constant<GJKType, GJKType::SIGNEDVOLUME>;
    using SP = std::integral_constant<GJKType, GJKType::SPHERE>;
    if(npType == NarrowPhaseType::SPHERE)
        std::forward<Func>(f)(SP{}, std::false_type{});
    else if(npType == NarrowPhaseType::GJK && !gjkAcc)
        std::forward<Func>(f)(J{}, std::false_type{});
    else if(npType == NarrowPhaseType::GJK && gjkAcc)
        std::forward<Func>(f)(J{}, std::true_type{});
    else if(npType == NarrowPhaseType::GJK_SV && !gjkAcc)
        std::forward<Func>(f)(SV{}, std::false_type{});
    else
        std::forward<Func>(f)(SV{}, std::true_type{});
}

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
    : m_neighborList(NeighborListFactory<T, M>::create(
          rigidBody, positions, orientations, CD, nObstacles, nParticles))
{
    size_t freeMem;
    if constexpr(M == MemType::HOST)
        freeMem = getAvailableHostMemory();
    else
        freeMem = getAvailableDeviceMemory();

    const auto& LCD       = GrainsParameters<T>::m_collisionDetection.linkedCellParameters;
    size_t estimatedPairs = static_cast<size_t>(nParticles) * LCD.initialNumberOfPairsPerParticle;
    size_t maxPairs       = static_cast<size_t>(nObstacles) * nParticles
                      + static_cast<size_t>(nParticles) * (nParticles - 1) / 2;
    estimatedPairs = std::min(estimatedPairs, maxPairs);

    m_pairBufferSize   = estimatedPairs;
    size_t sizePerPair = sizeof(uint8_t) + sizeof(uint);
    size_t sizeNeeded  = estimatedPairs * sizePerPair;
    GAssert(sizeNeeded < freeMem,
            "Not enough memory for pair-dependent buffers in CollisionDetectionModule!");

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
const uint2* CollisionDetectionModule<T, M>::getPairList() const
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
// For DEVICE + BV-ON, returns 0 so ComponentManager starts with zero-capacity buffers
// (grown lazily to bvPassCount on first run).
template <typename T, MemType M>
size_t CollisionDetectionModule<T, M>::getPairBufferSize() const
{
    if constexpr(M == MemType::DEVICE)
    {
        if(GrainsParameters<T>::m_collisionDetection.boundingVolumeType != BoundingVolumeType::OFF)
            return 0;
    }
    return m_pairBufferSize;
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
void CollisionDetectionModule<T, M>::run(GrainsMemBuffer<RigidBody<T>*, M>&  rigidBody,
                                         GrainsMemBuffer<Vector3<T>, M>&     positions,
                                         GrainsMemBuffer<Quaternion<T>, M>&  orientations,
                                         GrainsMemBuffer<Kinematics<T>, M>&  velocities,
                                         GrainsMemBuffer<Torce<T>, M>&       torces,
                                         GrainsMemBuffer<uint, M>&           bodyTags,
                                         GrainsMemBuffer<Vector3<T>, M>&     localPos,
                                         GrainsMemBuffer<Quaternion<T>, M>&  localQuat,
                                         GrainsMemBuffer<uint, M>&           masterSlot,
                                         GrainsMemBuffer<ContactInfo<T>, M>& contactInfo,
                                         GrainsMemBuffer<uint2, M>&          pairList,
                                         ComponentCounts&                    counts)
{
    auto& gt = GrainsParameters<T>::m_cdmTimer;
    // Flush any async GPU work from the previous iteration (e.g. MoveParticles launched after
    // the last CDM run) *before* CDM sub-stage timing starts.  This prevents that work from
    // being wrongly attributed to the NeighborList stage via the cudaStreamSynchronize(0)
    // inside prepareLinkedCellUpdate().  The sync is a no-op when timings are disabled.
    if constexpr(M == MemType::DEVICE)
        if(gt.isEnabled())
            cudaDeviceSynchronize();
    gt.start(CDMStage::Total);
    sortParticles(rigidBody,
                  positions,
                  orientations,
                  velocities,
                  torces,
                  bodyTags,
                  localPos,
                  localQuat,
                  masterSlot,
                  counts);
    updateNeighborList(positions, pairList, contactInfo, counts);
    detectCollisionsComponents(rigidBody.getData(),
                               positions,
                               orientations,
                               pairList,
                               contactInfo,
                               pairList,
                               bodyTags,
                               counts);
    gt.stop(CDMStage::Total);
}

// -------------------------------------------------------------------------------------------------
// Resizes per-pair buffers and ComponentManager's contactInfo / pairList
template <typename T, MemType M>
void CollisionDetectionModule<T, M>::resizePairBuffers(
    GrainsMemBuffer<uint2, M>& pairList, GrainsMemBuffer<ContactInfo<T>, M>& contactInfo, uint size)
{
    if constexpr(M == MemType::DEVICE)
    {
        m_bvPassFlags.grow(size);
        m_bvPassPairIndices.grow(size);
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
        m_cubTempStorageBytes = std::max(m_cubTempStorageBytes, newBytes);
        m_cubTempStorage.grow(m_cubTempStorageBytes);
        // For DEVICE + BV-ON, CM's contactInfo/pairList are grown lazily to bvPassCount after GJK
        const BoundingVolumeType bvType
            = GrainsParameters<T>::m_collisionDetection.boundingVolumeType;
        if(bvType == BoundingVolumeType::OFF)
        {
            contactInfo.grow(size);
            pairList.grow(size);
        }
    }
    else
    {
        contactInfo.grow(size);
        pairList.grow(size);
    }
}

// -------------------------------------------------------------------------------------------------
// Updates the neighbor list; resizes buffers and increments global counter when rebuilt
template <typename T, MemType M>
void CollisionDetectionModule<T, M>::updateNeighborList(
    GrainsMemBuffer<Vector3<T>, M>&     positions,
    GrainsMemBuffer<uint2, M>&          pairList,
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
        // For DEVICE + BV-ON: CDM reads NL directly in kernels; CM's pairList is grown lazily
        // after GJK to bvPassCount. For HOST or DEVICE + BV-OFF: copy NL to CM's pairList now.
        if constexpr(M == MemType::HOST)
        {
            m_neighborList->getBuffer().copyTo(pairList);
            numPairs = newSize;
        }
        else
        {
            const BoundingVolumeType bvType
                = GrainsParameters<T>::m_collisionDetection.boundingVolumeType;
            if(bvType == BoundingVolumeType::OFF)
            {
                m_neighborList->getBuffer().copyTo(pairList);
                numPairs = newSize;
            }
            // else: counts.numPairs will be set to bvPassCount after GJK
        }
        SS.neighborListUpdateCount++;
    }

    gt.stop(CDMStage::NeighborList);
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
    {
        ::buildShapeAndBVData(m_shapeData.getData(),
                              m_bvData.getData(),
                              rigidBodies,
                              repSlots.data(),
                              nUniqueShapes);
    }
}

// -------------------------------------------------------------------------------------------------
// Narrow-phase GJK using world-frame positions/quaternions directly.
// DEVICE BV-ON path writes compactly to ComponentManager's contactInfo and pairList buffers.
template <typename T, MemType M>
void CollisionDetectionModule<T, M>::detectCollisionsComponents(
    const RigidBody<T>* const*               rigidBodies,
    const GrainsMemBuffer<Vector3<T>, M>&    positions,
    const GrainsMemBuffer<Quaternion<T>, M>& orientations,
    const GrainsMemBuffer<uint2, M>&         pairList,
    GrainsMemBuffer<ContactInfo<T>, M>&      contactInfo,
    GrainsMemBuffer<uint2, M>&               cmPairList,
    const GrainsMemBuffer<uint, M>&          bodyTags,
    ComponentCounts&                         counts)
{
    const uint               nPairs  = m_neighborList->getSize();
    const BoundingVolumeType bvType  = GrainsParameters<T>::m_collisionDetection.boundingVolumeType;
    const NarrowPhaseType    npType  = GrainsParameters<T>::m_collisionDetection.narrowPhaseType;
    const bool               gjkAcc  = GrainsParameters<T>::m_collisionDetection.gjkAcceleration;
    const bool           usePrebuilt = GrainsParameters<T>::m_collisionDetection.usePrebuiltShapes;
    const uint*          bodyTagsData    = bodyTags.getData();
    const Vector3<T>*    positionData    = positions.getData();
    const Quaternion<T>* orientationData = orientations.getData();
    auto&                gt              = GrainsParameters<T>::m_cdmTimer;

    if constexpr(M == MemType::HOST)
    {
        const uint2* pairData = pairList.getData();
        // BVType is passed directly to detectCollisionsComponents_common, which forwards
        // it to closestPointsRigidBodies. The BV early-out and no-contact sentinel are handled
        // internally. Intra-composite pairs are skipped via an explicit guard before the call.
        auto run = [&](auto bvTypeTag) {
            constexpr BoundingVolumeType BVT = decltype(bvTypeTag)::value;
            gt.start(CDMStage::NarrowPhase);
            dispatchGJK(npType, gjkAcc, [&](auto gjkV_tag, auto gjkA_tag) {
                constexpr GJKType GJKV = decltype(gjkV_tag)::value;
                constexpr bool    GJKA = decltype(gjkA_tag)::value;

                auto detectPair = [&](uint pairID) {
                    if(usePrebuilt)
                        detectCollisionsComponents_common<T, GJKV, GJKA, BVT>(pairData,
                                                                              m_shapeData.getData(),
                                                                              bodyTagsData,
                                                                              positionData,
                                                                              orientationData,
                                                                              contactInfo.getData(),
                                                                              pairID);
                    else
                        detectCollisionsComponents_common<T, GJKV, GJKA, BVT>(pairData,
                                                                              rigidBodies,
                                                                              positionData,
                                                                              orientationData,
                                                                              contactInfo.getData(),
                                                                              pairID);
                };

                for(uint i = 0; i < nPairs; ++i)
                {
                    if(counts.numComposites > 0)
                    {
                        const uint2 p    = pairData[i];
                        const uint  tagA = bodyTagsData[p.x];
                        const uint  tagB = bodyTagsData[p.y];
                        if(isSubBody(tagA) && isSubBody(tagB)
                           && getCompositeIdx(tagA) == getCompositeIdx(tagB))
                        {
                            contactInfo.getData()[i].setOverlapDistance(T(1));
                            continue;
                        }
                    }
                    detectPair(i);
                }
            });
            gt.stop(CDMStage::NarrowPhase);
        };

        using OFFT = std::integral_constant<BoundingVolumeType, BoundingVolumeType::OFF>;
        using OBBT = std::integral_constant<BoundingVolumeType, BoundingVolumeType::OBB>;
        using OBCT = std::integral_constant<BoundingVolumeType, BoundingVolumeType::OBC>;
        if(bvType == BoundingVolumeType::OBB)
            run(OBBT{});
        else if(bvType == BoundingVolumeType::OBC)
            run(OBCT{});
        else
            run(OFFT{});
    }
    else
    {
        const uint2* nlPairData = m_neighborList->getBuffer().getData();

        auto launchDeviceGJK =
            [&](uint pairCountToProcess, const uint* activePairIndices, uint2* compactPairListOut) {
                if(pairCountToProcess == 0)
                    return;

                uint numBlocks, numThreads;
                computeOptimalThreadsAndBlocks(pairCountToProcess,
                                               GrainsParameters<T>::m_GPU,
                                               numBlocks,
                                               numThreads);

                dispatchGJK(npType, gjkAcc, [&](auto gjkV_tag, auto gjkA_tag) {
                    constexpr GJKType GJKV = decltype(gjkV_tag)::value;
                    constexpr bool    GJKA = decltype(gjkA_tag)::value;
                    if(usePrebuilt)
                        detectCollisionsComponents_Kernel<T, GJKV, GJKA>
                            <<<numBlocks, numThreads>>>(m_shapeData.getData(),
                                                        nlPairData,
                                                        bodyTagsData,
                                                        activePairIndices,
                                                        positionData,
                                                        orientationData,
                                                        contactInfo.getData(),
                                                        compactPairListOut,
                                                        pairCountToProcess);
                    else
                        detectCollisionsComponents_Kernel<T, GJKV, GJKA, BoundingVolumeType::OFF>
                            <<<numBlocks, numThreads>>>(rigidBodies,
                                                        nlPairData,
                                                        activePairIndices,
                                                        positionData,
                                                        orientationData,
                                                        contactInfo.getData(),
                                                        compactPairListOut,
                                                        pairCountToProcess);
                });
            };

        const bool needFilter = (bvType == BoundingVolumeType::OBB
                                 || bvType == BoundingVolumeType::OBC || counts.numComposites > 0);

        if(needFilter)
        {
            filterPairsBV(rigidBodies, positions, orientations, bodyTagsData, counts.numComposites);

            const uint bvPassCount = (uint)*m_bvPassPairCountMapped.getData();
            // Ensure output buffers are sized before writing compact world-frame GJK results.
            contactInfo.grow(bvPassCount);
            cmPairList.grow(bvPassCount);
            contactInfo.setSize(bvPassCount);
            cmPairList.setSize(bvPassCount);
            // GJK on compacted pairs (world frame, no relPos/relQuat needed)
            gt.start(CDMStage::NarrowPhase);
            launchDeviceGJK(bvPassCount, m_bvPassPairIndices.getData(), cmPairList.getData());
            gt.stop(CDMStage::NarrowPhase);
            counts.numPairs = bvPassCount;
        }
        else
        {
            // BV-OFF, no composites: run GJK over all pairs.
            gt.start(CDMStage::NarrowPhase);
            launchDeviceGJK(nPairs, nullptr, nullptr);
            gt.stop(CDMStage::NarrowPhase);
            counts.numPairs = nPairs;
        }
    }
}

// -------------------------------------------------------------------------------------------------
// BV filter using world-frame positions/quaternions; computes relative transforms on-the-fly
// inside filterPairsBV_Kernel to avoid a separate pre-pass. Reads NL pair list directly.
template <typename T, MemType M>
void CollisionDetectionModule<T, M>::filterPairsBV(
    const RigidBody<T>* const*               rigidBodies,
    const GrainsMemBuffer<Vector3<T>, M>&    positions,
    const GrainsMemBuffer<Quaternion<T>, M>& orientations,
    const uint*                              bodyTags,
    uint                                     numComposites)
{
    if constexpr(M == MemType::DEVICE)
    {
        auto& gt = GrainsParameters<T>::m_cdmTimer;
        gt.start(CDMStage::BVFilter);

        const uint   nPairs     = m_neighborList->getSize();
        const uint2* nlPairData = m_neighborList->getBuffer().getData();
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
        const bool        usePrebuilt = GrainsParameters<T>::m_collisionDetection.usePrebuiltShapes;
        const BVData<T>*  bvData      = m_bvData.getData();
        const Vector3<T>* positionData       = positions.getData();
        const Quaternion<T>* orientationData = orientations.getData();
        uint8_t*             bvPassFlags     = m_bvPassFlags.getData();

        auto launchBVFilter = [&](auto bvTypeTag) {
            constexpr BoundingVolumeType BVT = decltype(bvTypeTag)::value;
            if constexpr(BVT == BoundingVolumeType::OFF)
            {
                // BV is disabled: only run the OFF filter kernel when composites exist
                // (to cull intra-composite pairs). Otherwise, mark all pairs as pass.
                if(numComposites > 0)
                {
                    filterPairsBV_Kernel<T, BVT><<<numBlocks, numThreads>>>(rigidBodies,
                                                                            nlPairData,
                                                                            bodyTags,
                                                                            numComposites,
                                                                            positionData,
                                                                            orientationData,
                                                                            bvPassFlags,
                                                                            nPairs);
                }
                else
                {
                    cudaMemset(bvPassFlags, 1, static_cast<size_t>(nPairs) * sizeof(uint8_t));
                }
            }
            else if(usePrebuilt)
            {
                filterPairsBV_Kernel<T, BVT><<<numBlocks, numThreads>>>(bvData,
                                                                        nlPairData,
                                                                        bodyTags,
                                                                        numComposites,
                                                                        positionData,
                                                                        orientationData,
                                                                        bvPassFlags,
                                                                        nPairs);
            }
            else
            {
                filterPairsBV_Kernel<T, BVT><<<numBlocks, numThreads>>>(rigidBodies,
                                                                        nlPairData,
                                                                        bodyTags,
                                                                        numComposites,
                                                                        positionData,
                                                                        orientationData,
                                                                        bvPassFlags,
                                                                        nPairs);
            }
        };

        using OFFT = std::integral_constant<BoundingVolumeType, BoundingVolumeType::OFF>;
        using OBBT = std::integral_constant<BoundingVolumeType, BoundingVolumeType::OBB>;
        using OBCT = std::integral_constant<BoundingVolumeType, BoundingVolumeType::OBC>;

        if(bvType == BoundingVolumeType::OBB)
            launchBVFilter(OBBT{});
        else if(bvType == BoundingVolumeType::OBC)
            launchBVFilter(OBCT{});
        else
            launchBVFilter(OFFT{});

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
}
// -------------------------------------------------------------------------------------------------
// Sorts particles by Morton codes
template <typename T, MemType M>
void CollisionDetectionModule<T, M>::sortParticles(GrainsMemBuffer<RigidBody<T>*, M>& rigidBody,
                                                   GrainsMemBuffer<Vector3<T>, M>&    positions,
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

    // Contact-history tables are keyed by current slot indices. Reordering the arrays would
    // change slot ownership and make HookeMemory reuse tangential history for the wrong bodies.
    // Until contact history is keyed by a stable per-component ID, Morton reordering must stay
    // off whenever memory-enabled contact models are active.
    if(LC.sortFrequency > 0 && GP::m_isContactWithMemory)
    {
        SS.particlesSorted = false;
        return;
    }

    // Adaptive skin resizes and recenters the live linked-cell grid during neighbor-list updates.
    // ParticleSorter currently owns a separate Morton Cells object built once at construction time,
    // so its cell geometry drifts out of sync as soon as adaptive skin changes the linked-cell
    // size. Until the sorter is wired to the live linked-cell cell size, sorting must stay off
    // when adaptive skin is enabled.
    if(LC.sortFrequency > 0 && LC.updateFrequency > 0)
    {
        SS.particlesSorted = false;
        return;
    }

    if(LC.sortFrequency > 0 && SS.neighborListUpdateCount % LC.sortFrequency == 0)
    {
        if(m_particleSorter == nullptr)
            m_particleSorter = std::make_unique<ParticleSorter<T, M>>(nObstacles, nParticles);

        gt.start(CDMStage::Sort);
        m_particleSorter->sortParticles(rigidBody,
                                        positions,
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
                cudaErrCheck(cudaGetLastError());
                cudaErrCheck(cudaDeviceSynchronize());
            }
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
template class CollisionDetectionModule<float, MemType::DEVICE>;
template class CollisionDetectionModule<double, MemType::DEVICE>;
