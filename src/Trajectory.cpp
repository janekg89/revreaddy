/* Trajectory.cpp
 * 
 */

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

void Trajectory::record(std::vector<Particle> activeParticles, unsigned long int t)
{
	std::vector<std::array<double,3>>	currentCoordinates;
	for (auto&& particle : activeParticles)
	{
		currentCoordinates.push_back(particle.position);
	}	
	this->trajectory.push_back(currentCoordinates);
}

void Trajectory::writeBufferToFile()
{
	std::ofstream file;
	file.open("trajectory.ixyz", std::ofstream::out | std::ofstream::app);
	for (auto&& particles : this->trajectory)
	{
		file << particles.size() << "\n";
		file << "#comment" << "\n"; 
		for (auto&& coordinates : particles)
		{
			file << "type" << "\t";
			file << coordinates[0] << "\t";
			file << coordinates[1] << "\t";
			file << coordinates[2] << "\n";
		}
	}
	file.close();
}
