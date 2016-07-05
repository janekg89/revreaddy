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

TEST_F(FractionalTest, runfractional) {
    World w;
    Config c;
    FractionalDiffusion sim = FractionalDiffusion(&w, &c);
    c.boxsize = 10.;
    const unsigned long maxtime = 1000;
    sim.run(maxtime);


    // Random * r;

//std::cout <<  r  << std::endl;
    EXPECT_EQ(1, 1);


}
// have one scenario to do a pure reaction setup, a pure diffusion setup,
// an observable setup, a geometry setup