#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <Config.h>

namespace {
class ConfigUnittest : public ::testing::Test {};
}//namespace

TEST_F(ConfigUnittest, ctor) {
	Config c;
	EXPECT_EQ(c.timestep, 0.001);
	EXPECT_EQ(c.kT, 1.);
	EXPECT_EQ(c.isPeriodic, true);
	EXPECT_EQ(c.boxsize, 10.);
}

TEST_F(ConfigUnittest, new_Type) {
	Config c;
	// TODO check exceptions
	c.new_Type("superparticle", 2., 3.);
	EXPECT_EQ(c.particleTypes[0].name, "superparticle");
	EXPECT_EQ(c.particleTypes[0].radius, 2.);
	EXPECT_EQ(c.particleTypes[0].diffusionConstant, 3.);
}

TEST_F(ConfigUnittest, getNumberOfParticleTypes) {
	Config c;
	c.new_Type("testparticle", 12., 100.);
	c.new_Type("testparticle2", 12., 100.);
	EXPECT_EQ(c.getNumberOfParticleTypes(), 2);
}

TEST_F(ConfigUnittest, getParticleTypeStuff) {
	Config c;
	c.new_Type("ultraparticle", 1., 4.);
	EXPECT_EQ(c.getParticleTypeName(0), "ultraparticle");
	EXPECT_EQ(c.getParticleTypeRadius(0), 1.);
	EXPECT_EQ(c.getParticleTypeDiffusionConstant(0), 4.);
}

TEST_F(ConfigUnittest, new_Wall) {
	Config c;
	std::vector<double> normal = {1.,0.,0.};
	std::vector<double> point = {50.,0.,-10.};
	double strength = 5.;
	std::vector<unsigned> particleTypeIds = {0,2};
	c.new_Type("testparticle0", 0., 0.);
	c.new_Type("testparticle1", 0., 0.);
	c.new_Type("testparticle2", 0., 0.);
	c.new_Wall(normal, point, strength, particleTypeIds);
	EXPECT_EQ(c.geometries.size(), 1);
	EXPECT_THAT(c.geometries[0]->particleTypeIds, ::testing::ElementsAre(0,2));
	// cannot test for normal, point or strength since they are
	// Wall properties, but we have only access to Geometry properties
}

TEST_F(ConfigUnittest, new_DoubleWellZ) {
	Config c;
	double distanceMinima = 2.;
	double strength = 3.;
	std::vector<unsigned> particleTypeIds = {0,2};
	c.new_Type("testparticle0", 0., 0.);
	c.new_Type("testparticle1", 0., 0.);
	c.new_Type("testparticle2", 0., 0.);
	c.new_DoubleWellZ(distanceMinima, strength, particleTypeIds);
	EXPECT_EQ(c.geometries.size(), 1);
	EXPECT_THAT(c.geometries[0]->particleTypeIds, ::testing::ElementsAre(0,2));
	// cannot test for distanceMinima or strength since they are
	// DoubleWellZ properties
}

TEST_F(ConfigUnittest, deleteAllGeometries) {
	Config c;
	std::vector<unsigned> particleTypeIds = {0};
	std::vector<double> x = {1.,0.,0.};
	c.new_Type("testparticle", 1.,4.);
	c.new_DoubleWellZ(2., 4., particleTypeIds);
	c.new_Wall(x, x, 3., particleTypeIds);
	EXPECT_EQ(c.geometries.size(), 2);
	c.deleteAllGeometries();
	EXPECT_EQ(c.geometries.size(), 0);
}

TEST_F(ConfigUnittest, new_SoftRepulsion) {
	Config c;
	std::array<unsigned,2> affectedTuple = {0,0};
	double repulsionStrength = 2.;
	c.new_Type("testparticle", 2., 4.);
	c.new_SoftRepulsion("softrep", affectedTuple, repulsionStrength);
	EXPECT_EQ(c.interactions.size(), 1);
	EXPECT_EQ(c.interactions[0]->name, "softrep");
	EXPECT_EQ(c.interactions[0]->type, "SoftRepulsion");
	EXPECT_THAT(c.interactions[0]->affectedTuple, ::testing::ElementsAre(0,0));
	EXPECT_EQ(c.interactions[0]->cutoff, 4.);
	EXPECT_EQ(c.interactions[0]->parameters[0], 2.);
}

TEST_F(ConfigUnittest, new_LennardJones) {
	Config c;
	std::array<unsigned,2> affectedTuple = {0,0};
	double epsilon = 420.;
	c.new_Type("testparticle", 3., 4.);
	c.new_LennardJones("LJay", affectedTuple, epsilon);
	EXPECT_EQ(c.interactions.size(), 1);
	EXPECT_EQ(c.interactions[0]->name, "LJay");
	EXPECT_EQ(c.interactions[0]->type, "LennardJones");
	EXPECT_THAT(c.interactions[0]->affectedTuple, ::testing::ElementsAre(0,0));
	EXPECT_EQ(c.interactions[0]->cutoff, 2.5*(3.+3.));
	EXPECT_EQ(c.interactions[0]->parameters[0], 420.);
}

TEST_F(ConfigUnittest, deleteAllInteractions) {
	Config c;
	std::array<unsigned,2> affectedTuple = {0,0};
	double epsilon = 420.;
	c.new_Type("testparticle", 3., 4.);
	c.new_LennardJones("LJay", affectedTuple, epsilon);
	c.new_LennardJones("LJay", affectedTuple, epsilon);
	EXPECT_EQ(c.interactions.size(), 2);
	c.deleteAllInteractions();
	EXPECT_EQ(c.interactions.size(), 0);
}