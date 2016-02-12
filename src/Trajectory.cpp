/* Trajectory.cpp */

// TODO recPeriod and clearPeriod should govern when traj. is recorded and when
// it is written to file
#include "Trajectory.h"

Trajectory::Trajectory(unsigned long inRecPeriod, unsigned long inClearPeriod, std::string inFilename)
{
	this->recPeriod = inRecPeriod;
	this->clearPeriod = inClearPeriod;
	this->filename = inFilename;
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
		pt.particleCoordinates = particle.position;
		pt.particleTime = t;
		currentCoordinates.push_back(pt);
	}	
	this->trajectory.push_back(currentCoordinates);
}

void Trajectory::writeBufferToFile()
{
	std::ofstream file;
	file.open(this->filename, std::ofstream::out | std::ofstream::app);
	for (auto&& particles : this->trajectory) {
		file << particles.size() << "\n";
		file << "#timestep " << particles[0].particleTime << "\n"; 
		for (auto&& particle : particles) {
			file << "T" << particle.particleTypeId << "\t";
			file << particle.particleCoordinates[0] << "\t";
			file << particle.particleCoordinates[1] << "\t";
			file << particle.particleCoordinates[2] << "\n";
		}
	}
	file.close();
}