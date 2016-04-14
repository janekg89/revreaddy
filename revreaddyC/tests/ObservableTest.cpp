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

TEST_F(ObservableTest, EnergyRunAndSaveToFile) {
	// if files exists, delete them first
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

TEST_F(ObservableTest, AcceptanceRunAndSaveToFile) {
	// if files exists, delete them first
	boost::filesystem::path pathDiffDat("acceptance_diffusion_test.dat");
	boost::filesystem::path pathDiffH5("acceptance_diffusion_test.h5");
	boost::filesystem::path pathReactDat("acceptance_reactions_test.dat");
	boost::filesystem::path pathReactH5("acceptance_reactions_test.h5");
	boost::filesystem::remove(pathDiffH5);
	boost::filesystem::remove(pathDiffDat);
	boost::filesystem::remove(pathReactH5);
	boost::filesystem::remove(pathReactDat);
	// now actual test. with default impl, the acceptance is always = 1
	World w; Config c; Simulation s(&w, &c, ""); 
	s.new_Acceptance(1, "acceptance_reactions_test.dat", true);
	s.new_Acceptance(1, "acceptance_reactions_test.h5", true);
	s.new_Acceptance(1, "acceptance_diffusion_test.dat", false);
	s.new_Acceptance(1, "acceptance_diffusion_test.h5", false);
	// no particles added. just run for 4 timesteps -> there should be 5 entries in the observable
	s.run(4);
	s.writeAllObservablesToFile();
	// check h5 files
	hid_t reactFileId = H5Fopen("acceptance_reactions_test.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
	hid_t diffFileId = H5Fopen("acceptance_diffusion_test.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
	hsize_t reactDims[1], diffDims[1];
	H5LTget_dataset_info(reactFileId, "/acceptanceProbs", reactDims, NULL, NULL);
	H5LTget_dataset_info(diffFileId, "/acceptanceProbs", diffDims, NULL, NULL);
	EXPECT_EQ(reactDims[0], 5) << "expecting 5 entries in acceptance of reactions";
	EXPECT_EQ(diffDims[0], 5) << "expecting 5 entries in acceptance of diffusion";
	double data[5];
	H5LTread_dataset_double(reactFileId, "/acceptanceProbs", data);
	EXPECT_THAT(data, ::testing::ElementsAre(1,1,1,1,1));
	H5LTread_dataset_double(diffFileId, "/acceptanceProbs", data);
	EXPECT_THAT(data, ::testing::ElementsAre(1,1,1,1,1));
	H5Fclose(reactFileId);
	H5Fclose(diffFileId);
}

TEST_F(ObservableTest, MSDRunAndSaveToFile) {
	// if files exists, delete them first
	boost::filesystem::path pathH5("msd_test.h5");
	boost::filesystem::path pathDat("msd_test.dat");
	boost::filesystem::remove(pathH5);
	boost::filesystem::remove(pathDat);
	// now actual test. default implementation
	World w; Config c; Simulation s(&w, &c, "");
	c.new_Type("testparticle", 1., 1.);
	std::vector<double> pos = {0.,0.,0.};
	w.addParticle(pos, 0);
	s.new_MeanSquaredDisplacement(1, "msd_test.h5", 0);
	s.new_MeanSquaredDisplacement(1, "msd_test.dat", 0);
	s.run(4);
	s.writeAllObservablesToFile();
	// check h5 file
	hid_t fileId = H5Fopen("msd_test.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
	hsize_t dims[1];
	H5LTget_dataset_info(fileId, "/meanSquaredDisplacements", dims, NULL, NULL);
	EXPECT_EQ(dims[0], 5) << "expecting 5 entries in MSD";
	double msd[5];
	H5LTread_dataset_double(fileId, "/meanSquaredDisplacements", msd);
	EXPECT_EQ(msd[0], 0.) << "initially the displacement must be zero";
	EXPECT_NE(msd[4], 0.) << "the final value is usually non-zero";
	unsigned int numbers[5];
	H5LTread_dataset(fileId, "/numberOfParticles", H5T_NATIVE_UINT, numbers);
	EXPECT_THAT(numbers, ::testing::ElementsAre(1,1,1,1,1)) << "number of particles should be 1 for all times";
	H5Fclose(fileId);
}

TEST_F(ObservableTest, ParticleNumbersRunAndSaveToFile) {
	boost::filesystem::path pathH5("pNumbers_test.h5");
	boost::filesystem::path pathDat("pNumbers_test.dat");
	boost::filesystem::remove(pathH5);
	boost::filesystem::remove(pathDat);
	// actual test
	World w; Config c; Simulation s(&w, &c, "");
	c.new_Type("testparticle0", 1., 1.);
	c.new_Type("testparticle1", 1., 1.);
	std::vector<double> pos = {0.,0.,0.};
	w.addParticle(pos, 0);
	w.addParticle(pos, 0);
	w.addParticle(pos, 1);
	s.new_ParticleNumbers(1, "pNumbers_test.h5", 0);
	s.new_ParticleNumbers(1, "pNumbers_test.dat", 0);
	s.run(4);
	s.writeAllObservablesToFile();
	// check h5 file
	hid_t fileId = H5Fopen("pNumbers_test.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
	hsize_t dims[1];
	H5LTget_dataset_info(fileId, "/particleNumbers", dims, NULL, NULL);
	EXPECT_EQ(dims[0], 5) << "expecting 5 entries";
	unsigned long pNumbers[5];
	H5LTread_dataset(fileId, "/particleNumbers", H5T_NATIVE_ULONG, pNumbers);
	EXPECT_THAT(pNumbers, ::testing::ElementsAre(2,2,2,2,2)) << "number of 0-particles should be 2 for all times";
	H5Fclose(fileId);
}