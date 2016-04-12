/* Acceptance.cpp */

#include "Acceptance.h"

// TODO make this cleared automatically as well
Acceptance::Acceptance(unsigned int inRecPeriod, unsigned int inClearPeriod, std::string inFilename, bool inReactionsOrDiffusion) {
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

void Acceptance::record(World * world, double t) {
	if (this->reactionsOrDiffusion) {
		this->acceptanceProbs.push_back(world->acceptProbReactions);
	} else {
		this->acceptanceProbs.push_back(world->acceptProbDiffusion);		
	}
	this->times.push_back(t);
}

void Acceptance::writeToH5() {
	H5::H5File file(this->filename.c_str(), H5F_ACC_TRUNC);
	createExtendibleDataset(file, "acceptanceProbs", this->acceptanceProbs);
	createExtendibleDataset(file, "times", this->times);
	// upon writing clear the data buffer
	this->acceptanceProbs.clear();
	this->times.clear();
}

void Acceptance::writeToDat() {
	std::ofstream file;
	file.open(this->filename, std::ofstream::out);
	file << "Time\tAcceptanceProbability\n";
	for (unsigned int i=0; i<this->acceptanceProbs.size(); i++) {
		file << this->times[i] << "\t" << this->acceptanceProbs[i] << "\n";
	}
	file.close();
	// upon writing clear the data buffer
	this->acceptanceProbs.clear();
	this->times.clear();
}