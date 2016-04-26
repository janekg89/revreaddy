/* TrajectoryUnique.cpp */
#include "TrajectoryUnique.h"

TrajectoryUnique::TrajectoryUnique(unsigned long inRecPeriod, unsigned long inClearPeriod, std::string inFilename) {
	if (inClearPeriod < inRecPeriod) {
		throw Exception("The clearing-period of Trajectory must be larger or equal to the recording-period.");	
	}
	recPeriod = inRecPeriod;
	clearPeriod = inClearPeriod;
	clearedAutomatically = true;
	filename = inFilename;
	observableTypeName = "TrajectoryUnique";
	isSetup = false;
}

void TrajectoryUnique::setup(World * world, Config * config) {
	//TODO get uniqueIds
	numberOfParticles = world->particles.size();
	for (auto&& p : world->particles) {
		this->uniqueIds.push_back( p.uniqueId );
		this->typeIds.push_back( p.typeId );
	}
	for (auto i=0; i<uniqueIds.size(); ++i) {
		this->stillExists[i] = true;
	}
}

void TrajectoryUnique::record(World * world, double t) {
	Snapshot snapshot;
	snapshot.exists.resize(numberOfParticles);
	snapshot.time = t;
	snapshot.positions.resize(boost::extents[numberOfParticles][3]);
	snapshot.forces.resize(boost::extents[numberOfParticles][3]);
	for (auto i=0; i<uniqueIds.size(); ++i) {
		if ( stillExists[i] ) {
			auto index = this->findParticleIndex(world->particles, uniqueIds[i]);
			if (index == -1) {// particle not found
				stillExists[i] = false;
				snapshot.exists[i] = false;
			} else {// particle does exist
				snapshot.exists[i] = true;
				snapshot.positions[i][0] = world->particles[index].position[0];
				snapshot.positions[i][1] = world->particles[index].position[1];
				snapshot.positions[i][2] = world->particles[index].position[2];
				snapshot.forces[i][0] = world->particles[index].cumulativeForce[0];
				snapshot.forces[i][1] = world->particles[index].cumulativeForce[1];
				snapshot.forces[i][2] = world->particles[index].cumulativeForce[2];
			}
		} else {// particle is known to not exist
			snapshot.exists[i] = false;
		}
	}
	trajectory.push_back(snapshot);
}

void TrajectoryUnique::writeToDat() {
	this->writeToH5(); // no txt file implementation. create h5 instead
}

void TrajectoryUnique::writeToH5() {
	// if file exists, append to 'positions', 'forces', 'time' and 'exists' in first dimension
	// H5Fis_hdf5 returns positive value for true, and 0 for false
	boost::filesystem::path path( this->filename.c_str() );
	if ( boost::filesystem::exists(path) ) {
		this->appendToH5();
	}
	// if file does not exist, create it and write datasets to i
	// chunksize of file depends on clearPeriod. one chunk of position is [clearPeriod, N, 3]	 
	else {
		this->writeToNewH5();
	}
}

void TrajectoryUnique::appendToH5() {
	H5::H5File file(this->filename.c_str(), H5F_ACC_RDWR);
	auto T = trajectory.size();
	auto N = numberOfParticles;
	/* Buffer the time dependent data here */
	boost::multi_array<double,3> positions(boost::extents[T][N][3]);
	boost::multi_array<double,3> forces(boost::extents[T][N][3]);
	boost::multi_array<int,2> exists(boost::extents[T][N]);
	std::vector<double> times;
	bufferTimeDependentData(positions, forces, exists, times);
	/* Write the extended data */
	extendDataset(file, "positions", positions);
	extendDataset(file, "forces", forces);
	extendDataset(file, "exists", exists);
	extendDataset(file, "times", times);
	/* Finally clear the trajectory buffer */
	this->trajectory.clear();
}

void TrajectoryUnique::writeToNewH5() {
	H5::H5File file(this->filename.c_str(), H5F_ACC_TRUNC);
	auto T = trajectory.size();
	auto N = numberOfParticles;
	/* Buffer the time dependent data here. */
	boost::multi_array<double,3> positions(boost::extents[T][N][3]);
	boost::multi_array<double,3> forces(boost::extents[T][N][3]);
	boost::multi_array<int,2> exists(boost::extents[T][N]);
	std::vector<double> times;
	bufferTimeDependentData(positions, forces, exists, times);
	/* Write all time dependent datasets */
	createExtendibleDataset(file, "positions", positions);
	createExtendibleDataset(file, "forces", forces);
	createExtendibleDataset(file, "exists", exists);
	createExtendibleDataset(file, "times", times);
	/* Write time-independent information */
	createExtendibleDataset(file, "typeIds", this->typeIds);
	createExtendibleDataset(file, "uniqueIds", this->uniqueIds);
	/* Finally clear the trajectory buffer */
	this->trajectory.clear();
}

void TrajectoryUnique::bufferTimeDependentData(boost::multi_array<double,3>& positions, boost::multi_array<double,3>& forces, boost::multi_array<int,2>& exists, std::vector<double>& times) {
	for (auto i=0; i<trajectory.size(); ++i) {
		times.push_back( trajectory[i].time );
		for (auto j=0; j<numberOfParticles; ++j) {
			positions[i][j][0] = trajectory[i].positions[j][0];
			positions[i][j][1] = trajectory[i].positions[j][1];
			positions[i][j][2] = trajectory[i].positions[j][2];
			forces[i][j][0] = trajectory[i].forces[j][0];
			forces[i][j][1] = trajectory[i].forces[j][1];
			forces[i][j][2] = trajectory[i].forces[j][2];
			/* hdf5 does not support boolean. Therefore: 0 - false, 1 - true. */
			if (trajectory[i].exists[j]) {
				exists[i][j] = 1;
			} else {
				exists[i][j] = 0;
			}
		}
	}
}