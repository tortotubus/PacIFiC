#include "GrainsTimer.hh"

#include <gtest/gtest.h>

TEST(GrainsTimerTest, EnabledTimersRecordAndReset)
{
    GrainsSimTimer simulationTimer;
    simulationTimer.enable(true);
    simulationTimer.start(SimStage::Total);
    simulationTimer.stop(SimStage::Total);
    EXPECT_GE(simulationTimer[SimStage::Total], 0.0);
    simulationTimer.reset();
    EXPECT_DOUBLE_EQ(simulationTimer[SimStage::Total], 0.0);

    GrainsCDMTimer collisionTimer;
    collisionTimer.enable(true);
    collisionTimer.start(CDMStage::Total);
    collisionTimer.stop(CDMStage::Total);
    EXPECT_GE(collisionTimer[CDMStage::Total], 0.0);
    collisionTimer.reset();
    EXPECT_DOUBLE_EQ(collisionTimer[CDMStage::Total], 0.0);

    GrainsForceTimer forceTimer;
    forceTimer.enable(true);
    forceTimer.start(FMStage::Total);
    forceTimer.stop(FMStage::Total);
    EXPECT_GE(forceTimer[FMStage::Total], 0.0);
    forceTimer.reset();
    EXPECT_DOUBLE_EQ(forceTimer[FMStage::Total], 0.0);
}
