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
	// TODO check if still exists and record
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
	hsize_t T = trajectory.size();
	hsize_t N = numberOfParticles;
	// get the datasets
	H5::DataSet dsetPositions = file.openDataSet("positions");
	H5::DataSet dsetForces = file.openDataSet("forces");
	H5::DataSet dsetExists = file.openDataSet("exists");
	H5::DataSet dsetTimes = file.openDataSet("times");
	// get the filespaces
	H5::DataSpace fspacePositions(dsetPositions.getSpace());
	H5::DataSpace fspaceForces(dsetForces.getSpace());
	H5::DataSpace fspaceExists(dsetExists.getSpace());
	H5::DataSpace fspaceTimes(dsetTimes.getSpace());
	// old dims
	hsize_t oldPositionsForces[3];
	fspacePositions.getSimpleExtentDims(oldPositionsForces);
	hsize_t oldExists[2];
	fspaceExists.getSimpleExtentDims(oldExists);
	hsize_t oldTimes[1];
	fspaceTimes.getSimpleExtentDims(oldTimes);
	// extend the datasets by ext*
	hsize_t extPositionsForces[3] = {T, N, 3};
	hsize_t extExists[2] = {T, N};
	hsize_t extTimes[1] = {T};
	// new sizes are ... only enlarge first dimension T
	hsize_t newPositionsForces[3] = {
		oldPositionsForces[0] + extPositionsForces[0],
		oldPositionsForces[1],
		oldPositionsForces[2]
	};
	hsize_t newExists[2] = {
		oldExists[0] + extExists[0],
		oldExists[1]
	};
	hsize_t newTimes[1] = {oldTimes[0] + extTimes[0]};
	// extend the dsets to new dims and update the filespaces to new extent
	dsetPositions.extend(newPositionsForces);
	dsetForces.extend(newPositionsForces);
	dsetExists.extend(newExists);
	dsetTimes.extend(newTimes);
	fspacePositions = H5::DataSpace(dsetPositions.getSpace());
	fspaceForces = H5::DataSpace(dsetForces.getSpace());
	fspaceExists = H5::DataSpace(dsetExists.getSpace());
	fspaceTimes = H5::DataSpace(dsetTimes.getSpace());
	// select hyperslabs to be written to
	hsize_t offsetPositionsForces[3] = {oldPositionsForces[0], 0, 0};
	fspacePositions.selectHyperslab(H5S_SELECT_SET, extPositionsForces, offsetPositionsForces);
	fspaceForces.selectHyperslab(H5S_SELECT_SET, extPositionsForces, offsetPositionsForces);
	hsize_t offsetExists[2] = {oldExists[0], 0};
	fspaceExists.selectHyperslab(H5S_SELECT_SET, extExists, offsetExists);
	hsize_t offsetTimes[1] = {oldTimes[0]};
	fspaceTimes.selectHyperslab(H5S_SELECT_SET, extTimes, offsetTimes);
	// define memspaces from which extended data is read
	H5::DataSpace mspacePositions(3, extPositionsForces);
	H5::DataSpace mspaceForces(3, extPositionsForces, NULL);
	H5::DataSpace mspaceExists(2, extExists, NULL);
	H5::DataSpace mspaceTimes(1, extTimes, NULL);
	// buffer the time dependent data here
	boost::multi_array<double,3> positions(boost::extents[T][N][3]);
	boost::multi_array<double,3> forces(boost::extents[T][N][3]);
	boost::multi_array<bool,2> exists(boost::extents[T][N]);
	std::vector<double> times;
	bufferTimeDependentData(positions, forces, exists, times);
	// write the extended data
	dsetPositions.write(positions.data(), H5::PredType::NATIVE_DOUBLE, mspacePositions, fspacePositions);
	dsetForces.write(forces.data(), H5::PredType::NATIVE_DOUBLE, mspaceForces, fspaceForces);
	dsetExists.write(exists.data(), H5::PredType::NATIVE_HBOOL, mspaceExists, fspaceExists);
	dsetTimes.write(times.data(), H5::PredType::NATIVE_DOUBLE, mspaceTimes, fspaceTimes);

	// finally clear the trajectory buffer
	this->trajectory.clear();
}

void TrajectoryUnique::writeToNewH5() {
	H5::H5File file(this->filename.c_str(), H5F_ACC_TRUNC);
	hsize_t T = trajectory.size();
	hsize_t N = numberOfParticles;
	// define dimensions
	hsize_t dimsPositionsForces[3] = {T, N, 3}; // positions and forces are TxNx3, double
	hsize_t maxDimsPositionsForces[3] = {H5S_UNLIMITED, N, 3};
	hsize_t dimsExists[2] = {T, N}; // exists is TxN, bool
	hsize_t maxDimsExists[2] = {H5S_UNLIMITED, N};
	hsize_t dimsTimes[1] = {T}; // time is T, double
	hsize_t maxDimsTimes[1] = {H5S_UNLIMITED};
	hsize_t dimsIds[1] = {N}; // ids are N, unsigned and unsigned long long, these are only written once.
	// create dataspaces
	H5::DataSpace dspacePositionForces(3, dimsPositionsForces, maxDimsPositionsForces);
	H5::DataSpace dspaceExists(2, dimsExists, maxDimsExists);
	H5::DataSpace dspaceTimes(1, dimsTimes, maxDimsTimes);
	H5::DataSpace dspaceIds(1, dimsIds);
	// create datasets for positions and forces. therefore set properties to enable chunking.
	H5::DSetCreatPropList propPositionsForces;
	hsize_t chunkDimsPositionsForces[3] = {clearPeriod, N, 3};
	propPositionsForces.setChunk(3, chunkDimsPositionsForces);
	H5::DataSet dsetPositions = file.createDataSet("positions", H5::PredType::NATIVE_DOUBLE, dspacePositionForces, propPositionsForces);
	H5::DataSet dsetForces = file.createDataSet("forces", H5::PredType::NATIVE_DOUBLE, dspacePositionForces, propPositionsForces);
	// create datasets for time and exists. also enable chunking
	H5::DSetCreatPropList propTimes;
	hsize_t chunkDimsTimes[1] = {clearPeriod};
	propTimes.setChunk(1, chunkDimsTimes);
	H5::DSetCreatPropList propExists;
	hsize_t chunkDimsExists[2] = {clearPeriod, N};
	propExists.setChunk(2, chunkDimsExists);
	H5::DataSet dsetTimes = file.createDataSet("times", H5::PredType::NATIVE_DOUBLE, dspaceTimes, propTimes);
	H5::DataSet dsetExists = file.createDataSet("exists", H5::PredType::NATIVE_HBOOL, dspaceExists, propExists);
	// buffer the time dependent data here
	boost::multi_array<double,3> positions(boost::extents[T][N][3]);
	boost::multi_array<double,3> forces(boost::extents[T][N][3]);
	boost::multi_array<bool,2> exists(boost::extents[T][N]);
	std::vector<double> times;
	bufferTimeDependentData(positions, forces, exists, times);
	// write all time dependent datasets
	dsetPositions.write(positions.data(), H5::PredType::NATIVE_DOUBLE);
	dsetForces.write(forces.data(), H5::PredType::NATIVE_DOUBLE);
	dsetTimes.write(times.data(), H5::PredType::NATIVE_DOUBLE);
	dsetExists.write(exists.data(), H5::PredType::NATIVE_HBOOL);
	//create datasets for ids and write them. no chunking since these are fixed
	H5::DataSet dsetTypeIds = file.createDataSet("typeIds", H5::PredType::NATIVE_UINT, dspaceIds);
	H5::DataSet dsetUniqueIds = file.createDataSet("uniqueIds", H5::PredType::NATIVE_ULLONG, dspaceIds);
	dsetTypeIds.write(this->typeIds.data(), H5::PredType::NATIVE_UINT);
	dsetUniqueIds.write(this->uniqueIds.data(), H5::PredType::NATIVE_ULLONG);
	// finally clear the trajectory buffer
	this->trajectory.clear();
}

void TrajectoryUnique::bufferTimeDependentData(boost::multi_array<double,3>& positions, boost::multi_array<double,3>& forces, boost::multi_array<bool,2>& exists, std::vector<double>& times) {
	for (auto i=0; i<trajectory.size(); ++i) {
		times.push_back( trajectory[i].time );
		for (auto j=0; j<numberOfParticles; ++j) {
			positions[i][j][0] = trajectory[i].positions[j][0];
			positions[i][j][1] = trajectory[i].positions[j][1];
			positions[i][j][2] = trajectory[i].positions[j][2];
			forces[i][j][0] = trajectory[i].forces[j][0];
			forces[i][j][1] = trajectory[i].forces[j][1];
			forces[i][j][2] = trajectory[i].forces[j][2];
			exists[i][j] = trajectory[i].exists[j];
		}
	}
}