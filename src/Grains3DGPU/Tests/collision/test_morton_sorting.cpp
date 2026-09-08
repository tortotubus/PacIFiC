#include <gtest/gtest.h>

#include <array>

#include "BodyTag.hh"
#include "Box.hh"
#include "CollisionDetectionModule.hh"
#include "ContactInfo.hh"
#include "GrainsParameters.hh"
#include "Kinematics.hh"
#include "Quaternion.hh"
#include "RigidBody.hh"
#include "Torce.hh"
#include "Vector3.hh"

class MortonSortingHostTest : public ::testing::Test
{
protected:
    CollisionDetectionParameters<double> savedCollisionDetection;
    SimulationState<double>              savedSimulationState;
    bool                                 savedContactMemory = false;

    void SetUp() override
    {
        savedCollisionDetection = GrainsParameters<double>::m_collisionDetection;
        savedSimulationState    = GrainsParameters<double>::m_simulationState;
        savedContactMemory      = GrainsParameters<double>::m_isContactWithMemory;
    }

    void TearDown() override
    {
        GrainsParameters<double>::m_collisionDetection  = savedCollisionDetection;
        GrainsParameters<double>::m_simulationState     = savedSimulationState;
        GrainsParameters<double>::m_isContactWithMemory = savedContactMemory;
    }

    void configureGlobals(bool contactMemoryEnabled)
    {
        using GP = GrainsParameters<double>;

        GP::m_collisionDetection                    = CollisionDetectionParameters<double>{};
        GP::m_simulationState                       = SimulationState<double>{};
        GP::m_isContactWithMemory                   = contactMemoryEnabled;
        GP::m_collisionDetection.neighborListType   = NeighborListType::NSQ;
        GP::m_collisionDetection.boundingVolumeType = BoundingVolumeType::OFF;
        GP::m_collisionDetection.narrowPhaseType    = NarrowPhaseType::GJK;
        GP::m_collisionDetection.gjkAcceleration    = false;
        GP::m_collisionDetection.usePrebuiltShapes  = false;

        auto& LC                           = GP::m_collisionDetection.linkedCellParameters;
        LC.minCorner                       = Vector3<double>(0.0, 0.0, 0.0);
        LC.maxCorner                       = Vector3<double>(1.0, 1.0, 1.0);
        LC.minCellSize                     = 0.1;
        LC.cellSizeFactor                  = 1.0;
        LC.maxNumCellsPerObstacle          = 8;
        LC.initialNumberOfPairsPerParticle = 4;
        LC.updateFrequency                 = 0;
        LC.sortFrequency                   = 1;
    }

    std::array<double, 3> runCollisionDetection(bool contactMemoryEnabled)
    {
        configureGlobals(contactMemoryEnabled);

        auto* sharedBox = new Box<double>(0.02, 0.02, 0.02);

        GrainsMemBuffer<RigidBody<double>*, MemType::HOST> rigidBody;
        GrainsMemBuffer<Vector3<double>, MemType::HOST>    position;
        GrainsMemBuffer<Quaternion<double>, MemType::HOST> orientation;
        GrainsMemBuffer<Kinematics<double>, MemType::HOST> velocity;
        GrainsMemBuffer<Torce<double>, MemType::HOST>      torce;
        GrainsMemBuffer<uint, MemType::HOST>               bodyTag;
        GrainsMemBuffer<Vector3<double>, MemType::HOST>    localPos;
        GrainsMemBuffer<Quaternion<double>, MemType::HOST> localQuat;

        rigidBody.reserve(3);
        position.reserve(3);
        orientation.reserve(3);
        velocity.reserve(3);
        torce.reserve(3);
        bodyTag.reserve(3);
        localPos.reserve(3);
        localQuat.reserve(3);

        const Quaternion<double>    identity(0.0, 0.0, 0.0, 1.0);
        const std::array<double, 3> inputX = {0.75, 0.15, 0.45};

        for(uint i = 0; i < 3; ++i)
        {
            rigidBody.push_back(new RigidBody<double>(sharedBox, 0.001, 1000.0, 1));
            position.push_back(Vector3<double>(inputX[i], 0.0, 0.0));
            orientation.push_back(identity);
            velocity.push_back(Kinematics<double>());
            torce.push_back(Torce<double>());
            bodyTag.push_back(makeStandaloneBodyTag(0u));
            localPos.push_back(Vector3<double>(0.0, 0.0, 0.0));
            localQuat.push_back(identity);
        }

        GrainsMemBuffer<uint, MemType::HOST> masterSlot;
        masterSlot.initialize(1);
        masterSlot[0] = 0u;

        GrainsMemBuffer<ContactInfo<double>, MemType::HOST> contactInfo;
        GrainsMemBuffer<uint2, MemType::HOST>               pairList;
        contactInfo.initialize(1);
        pairList.initialize(1);

        ComponentCounts counts;
        counts.numObstacles = 0;
        counts.numParticles = 3;

        CollisionDetectionModule<double, MemType::HOST> cdm(
            &rigidBody,
            position,
            orientation,
            bodyTag,
            bodyTag.getData(),
            GrainsParameters<double>::m_collisionDetection,
            counts.numObstacles,
            counts.numParticles);

        cdm.run(rigidBody,
                position,
                orientation,
                velocity,
                torce,
                bodyTag,
                localPos,
                localQuat,
                masterSlot,
                contactInfo,
                pairList,
                counts);

        return {position[0][X], position[1][X], position[2][X]};
    }
};

TEST_F(MortonSortingHostTest, UsesLinkedCellCellSizeForMortonSorting)
{
    const auto sortedX = runCollisionDetection(false);

    EXPECT_DOUBLE_EQ(sortedX[0], 0.15);
    EXPECT_DOUBLE_EQ(sortedX[1], 0.45);
    EXPECT_DOUBLE_EQ(sortedX[2], 0.75);
    EXPECT_TRUE(GrainsParameters<double>::m_simulationState.particlesSorted);
}

TEST_F(MortonSortingHostTest, SkipsSortingWhenContactMemoryIsEnabled)
{
    const auto originalOrderX = runCollisionDetection(true);

    EXPECT_DOUBLE_EQ(originalOrderX[0], 0.75);
    EXPECT_DOUBLE_EQ(originalOrderX[1], 0.15);
    EXPECT_DOUBLE_EQ(originalOrderX[2], 0.45);
    EXPECT_FALSE(GrainsParameters<double>::m_simulationState.particlesSorted);
}