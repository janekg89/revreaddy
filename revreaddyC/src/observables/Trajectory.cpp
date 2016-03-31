/* Trajectory.cpp */

// TODO recPeriod and clearPeriod should govern when traj. is recorded and when
// it is written to file
#include "Trajectory.h"

Trajectory::Trajectory(unsigned long inRecPeriod, unsigned long inClearPeriod, std::string inFilename)
{
	if (inClearPeriod < inRecPeriod) {
		throw Exception("The clearing-period of Trajectory must be larger or equal to the recording-period.");
	}
	this->recPeriod = inRecPeriod;
	this->clearPeriod = inClearPeriod;
	this->filename = inFilename;
	this->nextWrittenTimeIncrement = 0;
}

Trajectory::~Trajectory() {}

/* Trajectory does not need configuration */
void Trajectory::configure(World * world, Config * config)
{
	return;
}

void Trajectory::record(World * world, double t)
{
	std::vector<particleTuple>	currentCoordinates;
	for (auto&& particle : world->particles) {
		particleTuple pt;
		pt.particleTypeId = particle.typeId;
		pt.particleUniqueId = particle.uniqueId;
		pt.particlePosition = particle.position;
		pt.particleForce = particle.cumulativeForce;
		pt.particleTime = t;
		currentCoordinates.push_back(pt);
	}	
	this->trajectory.push_back(currentCoordinates);
}

void Trajectory::writeBufferToH5() {
	// construct file and open it with "out" mode ("open for read and write")
	h5xx::file file(this->filename, h5xx::file::out);
	// TODO
}

void Trajectory::writeBufferToDat()
{
	std::ofstream file;
	file.open(this->filename, std::ofstream::out | std::ofstream::app);
	for (auto&& particles : this->trajectory) {
		file << particles.size() << "\n";
		file << "#timestep " << particles[0].particleTime << "\n"; 
		for (auto&& particle : particles) {
			file << "T" << particle.particleTypeId << "\t";
			file << particle.particlePosition[0] << "\t";
			file << particle.particlePosition[1] << "\t";
			file << particle.particlePosition[2] << "\n";
		}
	}
	file.close();
	// TODO clean trajectory
}