#include "GrainsPlatform.hh"
#include <vector>

#include "BodyTag.hh"
#include "CollisionDetectionCommon.hh"
#include "CollisionComputation.hh"
#include "CollisionComputationCPU.hh"
#include "CollisionComputationCUDA.hh"
#include "CollisionFiltering.hh"
#include "CollisionFilteringCUDA.hh"
#include "CollisionDetectionModule.hh"
#include "CollisionShapeData.hh"
#include "CollisionShapeDataCPU.hh"
#include "CollisionShapeDataCUDA.hh"
#include "ParticleSorting.hh"
#include "ParticleSortingCPU.hh"
#include "ParticleSortingCUDA.hh"
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
    if constexpr(M == MemType::DEVICE)
        m_filterBackend.state = grains::cuda::createCollisionFilterState<T>(estimatedPairs);

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
    if constexpr(M == MemType::DEVICE)
        grains::cuda::resizeCollisionFilterState(*m_filterBackend.state, size);

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
// Backends fill one entry per shape; the CUDA backend owns the temporary upload.
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

    const ShapeTableView<T> view{m_shapeData.getData(),
                                 m_bvData.getData(),
                                 rigidBodies,
                                 nUniqueShapes};
    if constexpr(M == MemType::HOST)
        grains::cpu::buildShapeTables(view, repSlots.data());
    else
        grains::cuda::buildShapeTables(view, repSlots.data());
}

// -------------------------------------------------------------------------------------------------
// Performs a narrow-phase GJK-based collision detection for each pair and writes contact info in
// A-local frame. DEVICE BV-ON and composite paths use the backend-owned selected pair list.
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
    if constexpr(M == MemType::DEVICE)
    {
        if(config.boundingVolumeType != BoundingVolumeType::OFF || counts.numComposites > 0)
        {
            const auto selected = filterPairsBV(rigidBodies,
                                                pairList,
                                                m_contactInfoLocal,
                                                bodyTags.getData(),
                                                counts.numComposites);
            count               = selected.count;
            active              = selected.indices;
        }
    }
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
    if constexpr(M == MemType::DEVICE)
    {
        if(config.boundingVolumeType != BoundingVolumeType::OFF || counts.numComposites > 0)
        {
            SelectedPairIndices selected;
            if(config.useRelativeTransformations)
            {
                if(config.boundingVolumeType != BoundingVolumeType::OFF)
                    computeRelativeTransformations(positions, orientations, pairList);
                selected = filterPairsBV(rigidBodies,
                                         pairList,
                                         contactInfo,
                                         bodyTags.getData(),
                                         counts.numComposites);
            }
            else
                selected = filterPairsBV_global(rigidBodies,
                                                positions,
                                                orientations,
                                                pairList,
                                                contactInfo,
                                                bodyTags.getData(),
                                                counts.numComposites);
            count  = selected.count;
            active = selected.indices;
        }
    }
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
// Delegates relative-frame BV/composite filtering and compaction to the CUDA backend.
template <typename T, MemType M>
SelectedPairIndices CollisionDetectionModule<T, M>::filterPairsBV(
    const RigidBody<T>* const*           rigidBodies,
    const GrainsMemBuffer<IndexPair, M>& pairList,
    GrainsMemBuffer<ContactInfo<T>, M>&  contactInfoLocal,
    const uint*                          bodyTags,
    uint                                 numComposites)
{
    if constexpr(M == MemType::DEVICE)
    {
        const auto&             config = GrainsParameters<T>::m_collisionDetection;
        auto&                   timer  = GrainsParameters<T>::m_cdmTimer;
        const BVPairView<T>     pairs{pairList.getData(),
                                  rigidBodies,
                                  m_bvData.getData(),
                                  bodyTags,
                                  contactInfoLocal.getData(),
                                  m_neighborList->getSize(),
                                  numComposites};
        const RelativeBVView<T> view{pairs, m_relPosition.getData(), m_relQuaternion.getData()};
        timer.start(CDMStage::BVFilter);
        const auto selected = grains::cuda::filterRelativePairs(
            *m_filterBackend.state,
            view,
            BVFilterOptions{config.boundingVolumeType, config.usePrebuiltShapes});
        timer.stop(CDMStage::BVFilter);
        return selected;
    }
    else
        return {};  // HOST retains its inline narrow-phase BV checks.
}

// -------------------------------------------------------------------------------------------------
// Delegates world-frame BV/composite filtering without a separate relative-transform pass.
template <typename T, MemType M>
SelectedPairIndices CollisionDetectionModule<T, M>::filterPairsBV_global(
    const RigidBody<T>* const*               rigidBodies,
    const GrainsMemBuffer<Vector3<T>, M>&    positions,
    const GrainsMemBuffer<Quaternion<T>, M>& orientations,
    const GrainsMemBuffer<IndexPair, M>&     pairList,
    GrainsMemBuffer<ContactInfo<T>, M>&      contactInfoLocal,
    const uint*                              bodyTags,
    uint                                     numComposites)
{
    if constexpr(M == MemType::DEVICE)
    {
        const auto&           config = GrainsParameters<T>::m_collisionDetection;
        auto&                 timer  = GrainsParameters<T>::m_cdmTimer;
        const BVPairView<T>   pairs{pairList.getData(),
                                  rigidBodies,
                                  m_bvData.getData(),
                                  bodyTags,
                                  contactInfoLocal.getData(),
                                  m_neighborList->getSize(),
                                  numComposites};
        const GlobalBVView<T> view{pairs, positions.getData(), orientations.getData()};
        timer.start(CDMStage::BVFilter);
        const auto selected = grains::cuda::filterGlobalPairs(
            *m_filterBackend.state,
            view,
            BVFilterOptions{config.boundingVolumeType, config.usePrebuiltShapes});
        timer.stop(CDMStage::BVFilter);
        return selected;
    }
    else
        return {};
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
            const uint           nTotal = nObstacles + nParticles;
            const MasterSlotView view{bodyTags.getData(), masterSlot.getData(), nTotal};
            if constexpr(M == MemType::HOST)
                grains::cpu::rebuildMasterSlots(view);
            else
                grains::cuda::rebuildMasterSlots<T>(view);
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
