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
// check that with k*tau = 1 a reaction is performed with same positions before and after
TEST_F(ReactionTest, EnzymaticCandidatesAndPerform) {
    // NOTE: s is now the Impl class
    World w; Config c; SimulationImpl s(&w, &c);
    c.new_Type("A", 1., 1.);
    c.new_Type("B", 1., 1.);
    c.new_Type("C", 1., 1.);
    // A +{2.} C {1.}<-->{1.} B +{2.} C
    c.new_Enzymatic("enz", 0, 1, 2, 1., 1., 2.);
    // position A and C with a distance of 1.5
    std::vector<double> posA = {1.,1.,1.};
    std::vector<double> posC = {-0.5,1.,1.};
    w.addParticle(posC, 2);
    w.addParticle(posA, 0);
    // calculate the distance (energy/forces) for the A-C pair. while doing so, a reactionEvent should
    // be created in world->reactionCandidates. Indices for both particles are 0 and 1.
    s.calculateSingleForceEnergyCheckReactionCandidate(0, 1);
    EXPECT_EQ(w.reactionCandidates.size(), 1) << "one reactionEvent should be present";
    EXPECT_EQ(w.reactionCandidates[0].forwardOrBackward, true) << "direction must be forward";
    EXPECT_THAT(w.reactionCandidates[0].participants, ::testing::ElementsAre(0, 1)) << "uniqueIds should be 0 and 1";
    EXPECT_EQ(w.reactionCandidates[0].reactionId, 0) << "reactionId is 0";
    // now perform reaction
    c.timestep = 1.;
    // tau = 1 and k = 1 --> reaction MUST be performed, since probability is approximated as p=k*tau
    s.propagateReactions();
    EXPECT_EQ(w.particles[1].typeId, 1) << "particle A must have changed to B";
    EXPECT_THAT(w.particles[1].position, ::testing::ElementsAre(1.,1.,1.)) << "the position mustn't change";
    EXPECT_EQ(w.particles[0].typeId, 2) << "catalyst particle still has type C";
    EXPECT_THAT(w.particles[0].position, ::testing::ElementsAre(-0.5,1.,1.)) << "catalyst particle is still at its "
                                                                                        "old pos";
}

TEST_F(ReactionTest, ConversionIsAffected) {
    std::vector<unsigned> forwardTypes = {0};
    std::vector<unsigned> backwardTypes = {1};
    Conversion conv("conv", forwardTypes, backwardTypes, 0.1, 0.2);
    EXPECT_TRUE(conv.isAffectedForward(0)) << "forward type is 0";
    EXPECT_TRUE(conv.isAffectedBackward(1)) << "backward type is 1";
    EXPECT_TRUE(conv.isAffectedForward(forwardTypes));
    EXPECT_TRUE(conv.isAffectedBackward(backwardTypes));
    EXPECT_FALSE(conv.isAffectedForward(backwardTypes));
    EXPECT_FALSE(conv.isAffectedBackward(forwardTypes));
    EXPECT_EQ(conv.forwardRate, 0.1);
    EXPECT_EQ(conv.backwardRate, 0.2);
}

TEST_F(ReactionTest, ConversionCandidatesAndPerform) {
    World w; Config c; SimulationImpl s(&w, &c);
    c.new_Type("A", 1., 1.);
    c.new_Type("B", 1., 1.);
    c.new_Conversion("convers", 0, 1, 1., 0.2);
    std::vector<double> posA = {0.,0.,0.};
    w.addParticle(posA, 0);
    s.setupUnimolecularCandidateTypes();
    EXPECT_EQ(s.unimolecularCandidateTypes.size(), 2);
    EXPECT_EQ(s.unimolecularCandidateTypes[0].reactionId, 0);
    EXPECT_EQ(s.unimolecularCandidateTypes[0].forwardOrBackward, true);
    EXPECT_EQ(s.unimolecularCandidateTypes[0].particleTypeId, 0);
    EXPECT_EQ(s.unimolecularCandidateTypes[1].reactionId, 0);
    EXPECT_EQ(s.unimolecularCandidateTypes[1].forwardOrBackward, false);
    EXPECT_EQ(s.unimolecularCandidateTypes[1].particleTypeId, 1);
    c.timestep = 1.;
    // tau = 1, k=1 --> conversion must be performed
    s.propagateReactions();
    EXPECT_THAT(w.particles[0].position, ::testing::ElementsAre(0.,0.,0.));
    EXPECT_EQ(w.particles[0].typeId, 1);
}

// the same for conversion and fusion