/* ParticleNumbers.cpp */

#include "ParticleNumbers.h"

ParticleNumbers::ParticleNumbers(
	unsigned long inRecPeriod,
	unsigned long inClearPeriod,
	std::string inFilename,
	unsigned inParticleTypeId)
{
	this->recPeriod = inRecPeriod;
	this->clearPeriod = inClearPeriod;
	this->filename = inFilename;
	this->particleTypeId = inParticleTypeId;
}

/* No configuration necessary */
void ParticleNumbers::configure(World * world, Config * config) {}

void ParticleNumbers::record(World * world,	double t)
{
	this->time.push_back(t);
	unsigned long counter = 0;
	for (unsigned i=0; i<world->particles.size(); i++) {
		if (world->particles[i].typeId == this->particleTypeId) {
			counter += 1;
		}
	}
	this->particleNumbers.push_back(counter);
}

void ParticleNumbers::writeBufferToFile()
{
	// first determine the file extension
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

void ParticleNumbers::writeBufferToH5()
{
	BinaryFile file(this->filename);
	file.addDatasetDouble("time", this->time);
	file.addDatasetUnsignedLong("particleNumbers", this->particleNumbers);
}

void ParticleNumbers::writeBufferToDat()
{
	std::ofstream file;
	file.open(this->filename);
	file << "Time\tParticleNumbers\n";
	for (unsigned i=0; i<particleNumbers.size(); i++) {
		file << this->time[i] << "\t";
		file << this->particleNumbers[i] << "\n";
	}
	file.close();
}