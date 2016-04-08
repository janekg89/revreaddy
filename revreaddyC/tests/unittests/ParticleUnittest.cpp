#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "Particle.h"
#include <vector>

namespace {
class ParticleUnittest : public ::testing::Test {};
}//namespace

TEST_F(ParticleUnittest, ctor) {
	Particle p;
	EXPECT_EQ(p.uniqueId, 0);
	EXPECT_EQ(p.typeId, 0);
	EXPECT_THAT(p.position, ::testing::ElementsAre(0.,0.,0.));
	EXPECT_THAT(p.boxCoordinates, ::testing::ElementsAre(0,0,0));
	EXPECT_THAT(p.cumulativeForce, ::testing::ElementsAre(0.,0.,0.));
	EXPECT_EQ(p.position.size(), p.position.capacity());
	EXPECT_EQ(p.boxCoordinates.size(), p.boxCoordinates.capacity());
	EXPECT_EQ(p.cumulativeForce.size(), p.cumulativeForce.capacity());
}

TEST_F(ParticleUnittest, move) {
	Particle p;
	std::vector<double> dx = {1.,1.,1.};
	p.position = {0.,0.,0.};
	p.move(dx);
	EXPECT_THAT(p.position, ::testing::ElementsAre(1.,1.,1.));
	dx = {-2.,+2.,0.};
	p.move(dx);
	EXPECT_THAT(p.position, ::testing::ElementsAre(-1.,3.,1.));
}

TEST_F(ParticleUnittest, addForce) {
	Particle p;
	std::vector<double> f = {1.,-2.,3.};
	p.addForce(f);
	EXPECT_THAT(p.cumulativeForce, ::testing::ElementsAre(1.,-2.,3.));
	f = {0.,5.,5000.};
	p.addForce(f);
	EXPECT_THAT(p.cumulativeForce, ::testing::ElementsAre(1.,3.,5003.));
}

TEST_F(ParticleUnittest, resetForce) {
	Particle p;
	p.cumulativeForce = {42.,0.,100.};
	p.resetForce();
	EXPECT_THAT(p.cumulativeForce, ::testing::ElementsAre(0.,0.,0.));
}