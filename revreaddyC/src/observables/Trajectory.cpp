/* Trajectory.cpp */

#include "Trajectory.h"

Trajectory::Trajectory(unsigned long inRecPeriod, unsigned long inClearPeriod, std::string inFilename)
{
	if (inClearPeriod < inRecPeriod) {
		throw Exception("The clearing-period of Trajectory must be larger or equal to the recording-period.");
	}
	this->clearedAutomatically = true;
	this->recPeriod = inRecPeriod;
	this->clearPeriod = inClearPeriod;
	this->filename = inFilename;
	this->observableTypeName = "Trajectory";
	this->nextRecordedTimeIncrement = 0;
	isSetup = false;
}

void Trajectory::record(World * world, double t) {
	Snapshot snapshot;
	snapshot.time = t;
	snapshot.timeIncrement = nextRecordedTimeIncrement;
	snapshot.numberOfParticles = world->particles.size();
	for (auto i = 0; i<world->particles.size(); ++i) {
		snapshot.positionsX.push_back( world->particles[i].position[0] );
		snapshot.positionsY.push_back( world->particles[i].position[1] );
		snapshot.positionsZ.push_back( world->particles[i].position[2] );
		snapshot.forcesX.push_back( world->particles[i].cumulativeForce[0] );
		snapshot.forcesY.push_back( world->particles[i].cumulativeForce[1] );
		snapshot.forcesZ.push_back( world->particles[i].cumulativeForce[2] );
		snapshot.typeIds.push_back( world->particles[i].typeId );
		snapshot.uniqueIds.push_back( world->particles[i].uniqueId );
	}
	trajectory.push_back(snapshot);
	nextRecordedTimeIncrement += 1;
}

void Trajectory::writeToH5() {
	this->writeToDat(); // currently no efficient h5 storage due to variable particle numbers
}

void Trajectory::writeToDat() {
	LOG_TRACE("Enter Trajectory::writeToDat")
	std::ofstream file;
	file.open(this->filename, std::ofstream::out | std::ofstream::app);
	for (auto&& snapshot : this->trajectory) {
		file << snapshot.numberOfParticles << "\n";
		file << "#time " << snapshot.time << "\n"; 
		for (auto i=0; i<snapshot.numberOfParticles; ++i) {
			file << "T" << snapshot.typeIds[i] << "\t";
			file << snapshot.positionsX[i] << "\t";
			file << snapshot.positionsY[i] << "\t";
			file << snapshot.positionsZ[i] << "\n";
		}
	}
	file.close();
	// empty trajectory after coordinates have been written
	this->trajectory.clear();
}