/* Acceptance.cpp */

#include "Acceptance.h"

Acceptance::Acceptance(
	unsigned int inRecPeriod,
	unsigned int inClearPeriod,
	std::string inFilename,
	bool inReactionsOrDynamics)
{
	this->recPeriod = inRecPeriod;
	this->clearPeriod = inClearPeriod;
	this->filename = inFilename;
	this->reactionsOrDynamics = inReactionsOrDynamics;
}

/* No configuration necessary */
void Acceptance::configure(World * world, Config * config) {}

void Acceptance::record(World * world, double t)
{
	if (this->reactionsOrDynamics) {
		this->acceptanceProbs.push_back(world->acceptProbReactions);
	}
	else {
		this->acceptanceProbs.push_back(world->acceptProbDynamics);		
	}

	this->times.push_back(t);
}

void Acceptance::writeBufferToFile()
{
	unsigned int lastDot = this->filename.find_last_of(".");
	std::string extension = this->filename.substr(lastDot);
	if ( (extension == ".h5") || (extension == ".hdf5") ) {
		this->writeBufferToH5();
	}
	else if ( (extension == ".dat") || (extension == ".txt") ) {
		this->writeBufferToDat();
	}
	else {
		this->writeBufferToDat();
	}
}

void Acceptance::writeBufferToH5()
{
	BinaryFile file(this->filename);
	file.addDatasetDouble(
		"time",
		this->times);
	file.addDatasetDouble(
		"acceptanceProb",
		this->acceptanceProbs);
}

void Acceptance::writeBufferToDat()
{
	std::ofstream file;
	file.open(this->filename, std::ofstream::out);
	file << "Time\tAcceptanceProbability\n";
	for (unsigned int i=0; i<this->acceptanceProbs.size(); i++) {
		file << this->times[i] << "\t" << this->acceptanceProbs[i] << "\n";
	}
	file.close();
}