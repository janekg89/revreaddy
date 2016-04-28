//
// Created by chris on 28.04.16.
//

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "Simulation.h"
#include "World.h"
#include "Config.h"

namespace  {
    class SimulationTest : public ::testing::Test{};
}//namespace

// TODO
// check that neighborlist is created and deleted correctly also when in between runs() something changes
TEST_F(SimulationTest, NeighborlistCreatedAndDeleted) {
    EXPECT_EQ(1,1);
}

// check correctly propagated periodic boundaries
// check distance calculation of Utils. fulfills minimum image convention
// correct forces and energies for soft, lj interaction of few particles
// check saving loading state/oldState
// check calculation of acceptances !!
// compare naive force calculation vs neighborlist calculation