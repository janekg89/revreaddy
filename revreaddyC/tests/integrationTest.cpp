//
// Created by chris on 23.05.16.
//

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "Simulation.h"

namespace {
    class integrationTest : public ::testing::Test{};
}

// have one scenario to do a pure reaction setup, a pure diffusion setup,
// an observable setup, a geometry setup