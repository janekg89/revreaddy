//
// Created by chris on 23.05.16.
//

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "Simulation.h"

namespace {
    class FractionalTest : public ::testing::Test {
    };
}

TEST_F(FractionalTest, singleParticleConversion) {
    World w;
    Config c;
    FractionalDiffusion sim(&w, &c);
    c.boxsize = 10.;
    c.new_Type("A", 1., 1.);
    c.new_Type("B", 1., 1.);
    c.new_Conversion("A->B", 0, 1, 1., 0.);
    c.timestep = 1.;
    w.alpha = 0.5;
    std::vector<double> pos = {0,0,0};
    w.addParticle(pos, 0);
    sim.run(5);
    EXPECT_EQ(w.particles.size(), 1);
    EXPECT_EQ(w.particles[0].typeId, 1);
    EXPECT_NE(w.particles[0].position[0], 0);
    EXPECT_NE(w.particles[0].position[1], 0);
    EXPECT_NE(w.particles[0].position[2], 0);
}
// have one scenario to do a pure reaction setup, a pure diffusion setup,
// an observable setup, a geometry setup