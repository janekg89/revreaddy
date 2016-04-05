/* ParticleNumbers.cpp */

#include "ParticleNumbers.h"

ParticleNumbers::ParticleNumbers(
	unsigned long inRecPeriod,
	unsigned long inClearPeriod,
	std::string inFilename,
	unsigned inParticleTypeId)
{
	this->recPeriod = inRecPeriod;
	this->clearedAutomatically = false;
	this->clearPeriod = inClearPeriod;
	this->filename = inFilename;
	observableTypeName = "ParticleNumbers";
	this->particleTypeId = inParticleTypeId;
	isSetup = false;
}

void ParticleNumbers::record(World * world,	double t) {
	this->time.push_back(t);
	unsigned long counter = 0;
	for (unsigned i=0; i<world->particles.size(); i++) {
		if (world->particles[i].typeId == this->particleTypeId) {
			counter += 1;
		}
	}
	this->particleNumbers.push_back(counter);
}

void ParticleNumbers::writeToH5() {
	BinaryFile file(this->filename);
	file.addDatasetDouble("time", this->time);
	file.addDatasetUnsignedLong("particleNumbers", this->particleNumbers);
}

void ParticleNumbers::writeToDat() {
	std::ofstream file;
	file.open(this->filename);
	file << "Time\tParticleNumbers\n";
	for (unsigned i=0; i<particleNumbers.size(); i++) {
		file << this->time[i] << "\t";
		file << this->particleNumbers[i] << "\n";
	}
	file.close();
}