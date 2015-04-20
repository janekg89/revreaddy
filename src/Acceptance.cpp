/* Acceptance.cpp */

#include "Acceptance.h"

Acceptance::Acceptance(
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

void Acceptance::record(
	std::vector<Particle>& activeParticles,
	double t)
{
	this->acceptanceProbs.push_back(this->simulation->currentAcceptance);
	this->times.push_back(t);
}

void Acceptance::writeBufferToFile()
{
	std::ofstream file;
	file.open(this->filename, std::ofstream::out);
	file << "Time\tAcceptanceProbability\n";
	for (unsigned int i=0; i<this->acceptanceProbs.size(); i++) {
		file << this->times[i] << "\t" << this->acceptanceProbs[i] << "\n";
	}
	file.close();
}
