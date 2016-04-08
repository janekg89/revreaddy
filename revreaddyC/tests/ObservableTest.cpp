#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "Simulation.h"
#include "World.h"
#include "Config.h"
#include <string>
#include <vector>
#include <iostream>
#include <boost/filesystem.hpp>

namespace {
class ObservableTest : public ::testing::Test{};
}//namespace

TEST_F(ObservableTest, TrajectoryUniqueRunAndSave) {
	// if file exists, delete it first
	boost::filesystem::path path("trajectoryunique_test.h5");
	bool removed = boost::filesystem::remove(path);
	if (removed) { std::cout << "testfile removed" << std::endl; }
	// now actual test
	World w;
	Config c;
	Simulation s(&w, &c, ""); // default impl
	c.new_Type("testparticle", 1., 1.);
	std::vector<double> pos = {0.,0.,0.};
	for (auto i=0; i<30; ++i) {
		w.addParticle(pos, 0);
	}
	s.new_TrajectoryUnique(1, 50, "trajectoryunique_test.h5");
	s.run(1000);
	s.writeAllObservablesToFile();
}