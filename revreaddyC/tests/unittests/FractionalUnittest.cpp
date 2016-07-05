//
// Created by janekg89 on 28.06.16.
//

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <Random.h>
#include <simulationImpls/generateIncrements.h>
#include <simulationImpls/FractionalDiffusion.h>

namespace {
    class FractionalUnittest : public ::testing::Test {};
}//namespace

TEST_F(FractionalUnittest, correctDimensions) {
    Random r;
    Random* rpoint= &r;
    //double ran_norm = r->normal();
    unsigned long N =1500;
    double D = 2.0;
    double tau = 0.5;
    double alpha = 0.7;
    //std::cout << rpoint->normal() << std::endl;
    auto result = janek::generateIncrements(N,  D, tau, alpha, rpoint );
    //std::cout <<  result.shape()[0] << "  "<< result.shape()[1]<<std::endl;
    EXPECT_EQ(result.shape()[0],  3);
    EXPECT_EQ(result.shape()[1],  N);

}

