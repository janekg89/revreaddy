#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "Random.h"

namespace {
class RandomUnittest : public ::testing::Test {};
}//namespace

TEST_F(RandomUnittest, ctor) {
	Random r1("mt19937");
	EXPECT_EQ(r1.getType(), "mt19937");
	Random r2("ranlxs0");
	EXPECT_EQ(r2.getType(), "ranlxs0");
	Random r3("taus");
	EXPECT_EQ(r3.getType(), "taus");
}

TEST_F(RandomUnittest, getNewSeed) {
	Random r;
	unsigned long int seed1 = r.getNewSeed();
	EXPECT_NE(seed1, r.getNewSeed()) << "two calls to getNewSeed should return different values";
}

TEST_F(RandomUnittest, normal) {
	Random r;
	double n2 = r.normal();
	double n1 = 0.;
	for (unsigned i=0;i<100;++i) {
		n1 = n2;
		n2 = r.normal();		
		EXPECT_NE(n1, n2) << "two subsequent random numbers should be different";
	}
}

TEST_F(RandomUnittest, normal3D) {
	Random r;
	std::vector<double> n1 = r.normal3D();
	EXPECT_EQ(n1.size(), 3);
	std::vector<double> n2 = {0.,0.,0.};
	for (unsigned i=0;i<100;++i) {
		n1 = n2;
		n2 = r.normal3D();		
		EXPECT_NE(n1[0], n2[0]) << "two subsequent random vectors should be different";
		EXPECT_NE(n1[1], n2[1]) << "two subsequent random vectors should be different";
		EXPECT_NE(n1[2], n2[2]) << "two subsequent random vectors should be different";
	}
}

TEST_F(RandomUnittest, uniform) {
	Random r;
	double u1 = r.uniform();
	EXPECT_GE(u1, 0.) << "uniform must be larger or equal zero";
	EXPECT_LT(u1, 1.) << "uniform must be smaller than one";
	double u2 = 0.;
	for (unsigned i=0;i<100;++i) {
		u1 = u2;
		u2 = r.uniform();		
		EXPECT_NE(u1, u2) << "two subsequent random numbers should be different";
		EXPECT_GE(u2, 0.) << "uniform must be larger or equal zero";
		EXPECT_LT(u2, 1.) << "uniform must be smaller than one";
	}
}