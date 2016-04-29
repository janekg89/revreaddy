//
// Created by chris on 28.04.16.
//

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "Simulation.h"
#include "utils.h"
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
TEST_F(SimulationTest, periodicBoundaryCondition) {
    World w; Config c; SimulationImpl s(&w, &c);
    c.new_Type("A", 1., 0.); // NOTE: zero diffusion coefficient
    c.boxsize = 10.;
    // place particle outside of the box
    std::vector<double> pos = {6., 0., 0.};
    w.addParticle(pos, 0);
    s.propagateDiffusion();
    EXPECT_THAT(w.particles[0].position, ::testing::ElementsAre(-4.,0.,0.)) << "the particle must now be at the other"
                                                                                       " end of the box";
}

/* Place two particles central in the box so that minimum distance is euclidean distance.
 * Also place two particles on the edges of the box so that minimum distance is the
 * distance "through the wall". */
TEST_F(SimulationTest, minimumImageConvention) {
    World w; Config c; SimulationImpl s(&w, &c);
    c.isPeriodic = true;
    c.boxsize = 10.;
    std::vector<double> x1 = {1.,1.,1.};
    std::vector<double> x2 = {-1.,-1.,-1.};
    std::vector<double> r = {0.,0.,0.}; // output value
    s.utils->getMinDistanceVector(r, x1, x2, c.isPeriodic, c.boxsize);
    EXPECT_THAT(r, ::testing::ElementsAre(-2.,-2.,-2.)) << "min distance should be euclidean distance";
    x1  = {4.,4.,4.};
    x2  = {-4.,-4.,-4.};
    s.utils->getMinDistanceVector(r, x1, x2, c.isPeriodic, c.boxsize);
    EXPECT_THAT(r, ::testing::ElementsAre(2.,2.,2.)) << "min distance should be distance 'through the wall'";
}

// compare naive force calculation vs neighborlist calculation
TEST_F(SimulationTest, naiveVsNeighborlistForces) {
    // create two simulations, worlds and configs. One uses neighborlist the other doesnt
    World w1; Config c1; SimulationImpl s1(&w1, &c1);
    c1.boxsize = 10.;
    c1.new_Type("A", 1.1, 1.);
    c1.new_SoftRepulsion("A<->A", {0,0}, 10.);
    World w2; Config c2; SimulationImpl s2(&w2, &c2);
    c2.boxsize = 10.;
    c2.new_Type("A", 1.1, 1.);
    c2.new_SoftRepulsion("A<->A", {0,0}, 10.);
    std::vector<double> pos = {0.,0.,0.};
    Random* r = s1.random;
    // place 64 particles at random positions, but identical in world1 and world2
    for (auto i=0; i<64 ; ++i) {
        pos = { r->uniform()*10.-5., r->uniform()*10.-5., r->uniform()*10.-5. };
        w1.addParticle(pos, 0);
        w2.addParticle(pos, 0);
    }
    s1.calculateInteractionForcesEnergiesNaive();
    s2.configureNeighborlist();
    s2.calculateInteractionForcesEnergiesWithLattice();
    // compare all forces of particles and total system energy
    for (auto i=0; i<64; ++i) {
        EXPECT_NEAR(w1.particles[i].cumulativeForce[0], w2.particles[i].cumulativeForce[0], 1e-11);
        EXPECT_NEAR(w1.particles[i].cumulativeForce[1], w2.particles[i].cumulativeForce[1], 1e-11);
        EXPECT_NEAR(w1.particles[i].cumulativeForce[2], w2.particles[i].cumulativeForce[2], 1e-11);
    }
    EXPECT_NEAR(w1.energy, w2.energy, 1e-11) << "the system energies must be the same";
    delete s2.neighborlist;
}

// correct forces and energies for soft, lj interaction of few particles
// check saving loading state/oldState
// check calculation of acceptances !!