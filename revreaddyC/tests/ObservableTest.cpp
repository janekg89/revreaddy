#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "Simulation.h"
#include "World.h"
#include "Config.h"
#include <string>
#include <vector>
#include <array>
#include <iostream>
#include <boost/filesystem.hpp>
#include <hdf5.h>
#include <hdf5_hl.h>

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

TEST_F(ObservableTest, EnergyRunAndSave) {
	// if file exists, delete it first
	boost::filesystem::path pathH5("energy_test.h5");
	boost::filesystem::path pathDat("energy_test.dat");
	boost::filesystem::remove(pathH5);
	boost::filesystem::remove(pathDat);
	// now actual test
	World w; Config c; Simulation s(&w, &c, ""); // default impl
	c.new_Type("testparticle", 0.14, 1.);
	std::array<unsigned,2> affected = {0, 0};
	c.new_SoftRepulsion("testrepulsion", affected, 0.1);
	std::vector<double> pos = {0.,0.,0.};
	for (auto i=0; i<30; ++i) {
		pos[0] = 0.1*i; // particles are lined up from 0. to 3. with distance 0.1
		w.addParticle(pos, 0);
	}
	s.new_Energy(1, "energy_test.h5");
	s.new_Energy(1, "energy_test.dat");
	s.run(100);
	s.writeAllObservablesToFile();
	hid_t file_id = H5Fopen("energy_test.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
	hsize_t dims[1];
	H5LTget_dataset_info(file_id, "/energies", dims, NULL, NULL);
	EXPECT_EQ(dims[0], 101) << "expecting the dataset to have 101 entries.";
	double data[101];
	H5LTread_dataset_double(file_id, "/energies", data);
	// initial energy value should be 0.11188
	EXPECT_NEAR(data[0], 0.11188, 1e-10) << "first energy value should be 0.11188";
	H5Fclose(file_id);
}