#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "Config.h"

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

TEST_F(ConfigUnittest, getNumberParticleTypes) {
	Config c;
	c.new_Type("testparticle", 12., 100.);
	c.new_Type("testparticle2", 12., 100.);
	EXPECT_EQ(c.getNumberParticleTypes(), 2);
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
	c.new_Wall("testwand",normal, point, strength, particleTypeIds);
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
	c.new_DoubleWellZ("doppelbrunnen", distanceMinima, strength, particleTypeIds);
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
	c.new_DoubleWellZ("doppelbrunnen", 2., 4., particleTypeIds);
	c.new_Wall("der wand", x, x, 3., particleTypeIds);
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

TEST_F(ConfigUnittest, new_Conversion) {
	Config c;
    c.new_Type("A", 1., 1.);
    c.new_Type("B", 1., 1.);
    // A <--> B
	c.new_Conversion("conv", 0, 1, 0.2, 0.5);
    EXPECT_EQ(c.reactions.size(), 1);
    EXPECT_THAT(c.reactions[0]->forwardTypes, ::testing::ElementsAre(0));
    EXPECT_THAT(c.reactions[0]->backwardTypes, ::testing::ElementsAre(1));
    EXPECT_EQ(c.reactions[0]->forwardRate, 0.2);
    EXPECT_EQ(c.reactions[0]->backwardRate, 0.5);
    EXPECT_EQ(c.reactions[0]->name, "conv");
    EXPECT_EQ(c.reactions[0]->type, "Conversion");
}

TEST_F(ConfigUnittest, new_Fusion) {
    Config c;
    c.new_Type("A", 1., 1.);
    c.new_Type("B", 1., 1.);
    c.new_Type("C", 1., 1.);
    // A + B <--> C
    c.new_Fusion("fusi", 0, 1, 2, 0.3, 0.1, 3.);
    EXPECT_EQ(c.reactions.size(), 1);
    EXPECT_THAT(c.reactions[0]->forwardTypes, ::testing::ElementsAre(0,1));
    EXPECT_THAT(c.reactions[0]->backwardTypes, ::testing::ElementsAre(2));
    EXPECT_EQ(c.reactions[0]->forwardRate, 0.3);
    EXPECT_EQ(c.reactions[0]->backwardRate, 0.1);
    EXPECT_EQ(c.reactions[0]->reactionDistance, 3.);
    EXPECT_EQ(c.reactions[0]->name, "fusi");
    EXPECT_EQ(c.reactions[0]->type, "Fusion");
}

TEST_F(ConfigUnittest, new_Enzymatic) {
    Config c;
    c.new_Type("A", 1., 1.);
    c.new_Type("B", 1., 1.);
    c.new_Type("C", 1., 1.);
    // A + C <--> B + C
    c.new_Enzymatic("enz", 0, 1, 2, 0.2, 0.6, 4.);
    EXPECT_EQ(c.reactions.size(), 1);
    EXPECT_THAT(c.reactions[0]->forwardTypes, ::testing::ElementsAre(0,2));
    EXPECT_THAT(c.reactions[0]->backwardTypes, ::testing::ElementsAre(1,2));
    EXPECT_EQ(c.reactions[0]->forwardRate, 0.2);
    EXPECT_EQ(c.reactions[0]->backwardRate, 0.6);
    EXPECT_EQ(c.reactions[0]->reactionDistance, 4.);
    EXPECT_EQ(c.reactions[0]->name, "enz");
    EXPECT_EQ(c.reactions[0]->type, "Enzymatic");
}

