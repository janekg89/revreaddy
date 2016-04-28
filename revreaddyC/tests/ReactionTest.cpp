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

// check isAffected forward and backward
TEST_F(ReactionTest, EnzymaticIsAffected) {
    std::vector<unsigned> forwardTypes = {0,2};
    std::vector<unsigned> backwardTypes = {1,2};
    double forwardRate = 0.2;
    double backwardRate = 0.3;
    double reactionDistance = 2.;
    Enzymatic enz("enz", forwardTypes, backwardTypes, forwardRate, backwardRate, reactionDistance);
    EXPECT_FALSE(enz.isAffectedBackward(0)) << "a single particle cannot be affected";
    EXPECT_FALSE(enz.isAffectedForward(1)) << "a single particle cannot be affected";
    EXPECT_TRUE(enz.isAffectedForward(forwardTypes));
    EXPECT_FALSE(enz.isAffectedForward(backwardTypes));
    EXPECT_TRUE(enz.isAffectedBackward(backwardTypes));
    EXPECT_FALSE(enz.isAffectedBackward(forwardTypes));
    EXPECT_EQ(enz.forwardRate, forwardRate);
    EXPECT_EQ(enz.backwardRate, backwardRate);
    EXPECT_EQ(enz.forwardTypes, forwardTypes);
    EXPECT_EQ(enz.backwardTypes, backwardTypes);
    EXPECT_EQ(enz.reactionDistance, reactionDistance);
}

// check that after calcSingle..(i,j) the pair appears in world->reactionCandidates. forward and backward
TEST_F(ReactionTest, EnzymaticCandidates) {
    // NOTE: s is now the Impl class
    World w; Config c; SimulationImpl s(&w, &c);
    c.new_Type("A", 1., 1.);
    c.new_Type("B", 1., 1.);
    c.new_Type("C", 1., 1.);
    // A +{2.} C {0.}<-->{0.} B +{2.} C
    c.new_Enzymatic("enz", 0, 1, 2, 0., 0., 2.);
    // position A and C with a distance of 1.5
    std::vector<double> posA = {1.,1.,1.};
    std::vector<double> posC = {-0.5,1.,1.};
    w.addParticle(posA, 0);
    w.addParticle(posC, 2);
    // calculate the distance (energy/forces) for the A-C pair. while doing so, a reactionEvent should
    // be created in world->reactionCandidates. Indices for both particles are 0 and 1.
    s.calculateSingleForceEnergyCheckReactionCandidate(0, 1);
    EXPECT_EQ(w.reactionCandidates.size(), 1) << "one reactionEvent should be present";
    EXPECT_EQ(w.reactionCandidates[0].forwardOrBackward, true) << "direction must be forward";
    EXPECT_THAT(w.reactionCandidates[0].participants, ::testing::ElementsAre(0, 1)) << "uniqueIds should be 0 and 1";
    EXPECT_EQ(w.reactionCandidates[0].reactionId, 0) << "reactionId is 0";
}

// TODO
// check that with k*tau = 1 a reaction is performed with same positions before and after
TEST_F(ReactionTest, EnzymaticPerform) {
    EXPECT_EQ(1,1);
}

// the same for conversion and fusion