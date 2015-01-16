/* Trajectory.cpp */

// TODO recPeriod and clearPeriod should govern when traj. is recorded and when
// it is written to file
#include "Trajectory.h"

Trajectory::Trajectory()
{
	this->recPeriod = 1;
	this->clearPeriod = 0;
}

Trajectory::~Trajectory()
{
	this->recPeriod = 0;
}

void Trajectory::record(
	std::vector<Particle> activeParticles,
	unsigned long int t)
{
	std::vector<particleTuple>	currentCoordinates;
	for (auto&& particle : activeParticles)
	{
		particleTuple pt;
		pt.particleType = particle.type;
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
	for (auto&& particles : this->trajectory)
	{
		file << particles.size() << "\n";
		file << "#timestep " << particles[0].particleTime << "\n"; 
		for (auto&& particle : particles)
		{
			file << particle.particleType << "\t";
			file << particle.particleCoordinates[0] << "\t";
			file << particle.particleCoordinates[1] << "\t";
			file << particle.particleCoordinates[2] << "\n";
		}
	}
	file.close();
}
