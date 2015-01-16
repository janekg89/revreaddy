/* Trajectory.cpp */
#include "TrajectorySingle.h"

void TrajectorySingle::record(
	std::vector<Particle> activeParticles,
	unsigned long int t)
{
	std::vector<double> currentCoordinates = {0.,0.,0.};
	Particle * particle = &activeParticles[0];
	currentCoordinates[0] = particle->position[0];
	currentCoordinates[1] = particle->position[1];
	currentCoordinates[2] = particle->position[2];
	this->trajectory.push_back(currentCoordinates);
}

void TrajectorySingle::writeBufferToFile()
{

}
