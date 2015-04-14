/* Energy.cpp */

#include "Energy.h"

Energy::Energy(
	unsigned int inRecPeriod,
	unsigned int inClearPeriod,
	Simulation * inSimulation,
	std::string inFilename)
{
	this->recPeriod = inRecPeriod;
	this->clearPeriod = inClearPeriod;
	this->simulation = inSimulation;
	this->filename = inFilename;
}

void Energy::record(
	std::vector<Particle>& activeParticles,
	double t)
{
	this->energies.push_back(this->simulation->energy);
	this->times.push_back(t);
}

void Energy::writeBufferToFile()
{
	std::ofstream file;
	file.open(this->filename, std::ofstream::out);
	file << "Time\tEnergy\n";
	for (unsigned int i=0; i<this->energies.size(); i++) {
		file << this->times[i] << "\t" << this->energies[i] << "\n";
	}
	file.close();
}
