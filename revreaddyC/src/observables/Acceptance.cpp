/* Acceptance.cpp */

#include "Acceptance.h"

Acceptance::Acceptance(
	unsigned int inRecPeriod,
	unsigned int inClearPeriod,
	std::string inFilename,
	bool inReactionsOrDiffusion)
{
	this->recPeriod = inRecPeriod;
	this->clearedAutomatically = false;
	this->clearPeriod = inClearPeriod;
	this->filename = inFilename;
	observableTypeName = "Acceptance";
	this->reactionsOrDiffusion = inReactionsOrDiffusion;
	isSetup = false;
}

/* No configuration necessary */
void Acceptance::configure(World * world, Config * config) {}

void Acceptance::record(World * world, double t)
{
	if (this->reactionsOrDiffusion) {
		this->acceptanceProbs.push_back(world->acceptProbReactions);
	}
	else {
		this->acceptanceProbs.push_back(world->acceptProbDiffusion);		
	}

	this->times.push_back(t);
}

void Acceptance::writeToH5()
{
	BinaryFile file(this->filename);
	file.addDatasetDouble(
		"time",
		this->times);
	file.addDatasetDouble(
		"acceptanceProb",
		this->acceptanceProbs);
}

void Acceptance::writeToDat()
{
	std::ofstream file;
	file.open(this->filename, std::ofstream::out);
	file << "Time\tAcceptanceProbability\n";
	for (unsigned int i=0; i<this->acceptanceProbs.size(); i++) {
		file << this->times[i] << "\t" << this->acceptanceProbs[i] << "\n";
	}
	file.close();
}