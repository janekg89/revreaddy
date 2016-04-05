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
}

void TrajectoryUnique::record(World * world, double t) {
	// TODO check if still exists and record
}

void TrajectoryUnique::writeToDat() {
	this->writeToH5();
}

void TrajectoryUnique::writeToH5() {
	// if file exists, append to 'positions', 'forces', 'time' and 'exists' in first dimension
	// if file does not exist, create it and write datasets to i
	// chunksize of file depends on clearPeriod. one chunk of position is [clearPeriod, N, 3]
}