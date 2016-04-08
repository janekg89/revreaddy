#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "World.h"
#include "Particle.h"

namespace {
class WorldUnittest : public ::testing::Test {};
}//namespace

TEST_F(WorldUnittest, ctor) {
	World w;
	EXPECT_EQ(w.cumulativeRuntime, 0.);
	EXPECT_EQ(w.energy, 0.);
	EXPECT_EQ(w.oldEnergy, 0.);
	EXPECT_EQ(w.acceptProbDiffusion, 1.);
	EXPECT_EQ(w.acceptProbReactions, 1.);
	EXPECT_EQ(w.acceptionsDiffusion, 0);
	EXPECT_EQ(w.rejectionsDiffusion, 0);
	EXPECT_EQ(w.acceptionsReactions, 0);
	EXPECT_EQ(w.rejectionsReactions, 0);
	EXPECT_EQ(w.uniqueIdCounter, 0);
}

TEST_F(WorldUnittest, addParticle) {
	World w;
	std::vector<double> x = {-2.,3.,1000.};
	unsigned id = 4;
	w.addParticle(x, id);
	EXPECT_THAT(w.particles[0].position, ::testing::ElementsAre(-2.,3.,1000.));
	EXPECT_EQ(w.particles[0].typeId, 4);
	EXPECT_EQ(w.particles[0].uniqueId, 0);
	EXPECT_EQ(w.uniqueIdCounter, 1);
	EXPECT_EQ(w.particles.size(), 1);
}

TEST_F(WorldUnittest, removeParticle) {
	World w;
	std::vector<double> x = {42.,1337.,-666.};
	unsigned id = 2000;
	w.addParticle(x, id);
	w.addParticle(x, id);
	EXPECT_EQ(w.particles.size(), 2);
	w.removeParticle(0);
	EXPECT_EQ(w.particles.size(), 1);
	w.removeParticle(0);
	EXPECT_EQ(w.particles.size(), 0);
	// TODO test if exception is thrown
}

TEST_F(WorldUnittest, getNumberOfParticles) {
	World w;
	std::vector<double> x = {42.,1337.,-666.};
	unsigned id = 2;
	w.addParticle(x, id);
	w.addParticle(x, id);
	w.addParticle(x, id);
	EXPECT_EQ(w.getNumberOfParticles(), 3);
	w.removeParticle(2);
	EXPECT_EQ(w.getNumberOfParticles(), 2);
	w.removeParticle(0);
	EXPECT_EQ(w.getNumberOfParticles(), 1);
}

TEST_F(WorldUnittest, getPosition) {
	World w;
	std::vector<double> x = {0.,2.,-50.};
	unsigned id = 2;
	w.addParticle(x, id);
	EXPECT_THAT(w.getPosition(0), ::testing::ElementsAre(0.,2.,-50.));
	// TODO test exception
}

TEST_F(WorldUnittest, setPosition) {
	World w;
	std::vector<double> x = {0.,2.,-50.};
	unsigned id = 5;
	w.addParticle(x, id);
	x = {5.,-10.,0.};
	w.setPosition(0, x);
	EXPECT_THAT(w.particles[0].position, ::testing::ElementsAre(5.,-10.,0.));
}

TEST_F(WorldUnittest, getTypeId) {
	World w;
	std::vector<double> x = {0.,2.,-50.};
	unsigned id = 5;
	unsigned id2 = 80;
	w.addParticle(x, id);
	w.addParticle(x, id2);
	EXPECT_EQ(w.getTypeId(0), 5);
	EXPECT_EQ(w.getTypeId(1), 80);
}

TEST_F(WorldUnittest, setTypeId) {
	World w;
	std::vector<double> x = {0.,2.,-50.};
	unsigned id = 5;
	unsigned id2 = 80;
	w.addParticle(x, id);
	w.addParticle(x, id2);
	w.setTypeId(0, 0);
	w.setTypeId(1, 200);
	EXPECT_EQ(w.particles[0].typeId, 0);
	EXPECT_EQ(w.particles[1].typeId, 200);
}

TEST_F(WorldUnittest, deleteAllParticles) {
	World w;
	std::vector<double> x = {0.,2.,-50.};
	unsigned id = 5;
	for (unsigned i=0; i<2000; ++i) {
		w.addParticle(x, id);
	}
	EXPECT_EQ(w.getNumberOfParticles(), 2000);
	w.deleteAllParticles();
	EXPECT_EQ(w.getNumberOfParticles(), 0);
}