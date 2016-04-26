//
// Created by chris on 22.04.16.
//

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "Simulation.h"
#include "World.h"
#include "Config.h"
#include "Enzymatic.h"

namespace {
    class ReactionTest : public ::testing::Test{};
}//namespace

TEST_F(ReactionTest, PerformForward) {
    // @todo
    EXPECT_EQ(1,1);
}
