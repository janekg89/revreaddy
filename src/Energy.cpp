/* Energy.cpp */

#include "Energy.h"

Energy::Energy(
	unsigned int inRecPeriod,
	unsigned int inClearPeriod,
	std::string inFilename)
{
	this->recPeriod = inRecPeriod;
	this->clearPeriod = inClearPeriod;
	this->filename = inFilename;
}

void Energy::record(
	World * world,
	double t)
{
	this->energies.push_back(world->energy);
	this->times.push_back(t);
}

void Energy::writeBufferToFile()
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

void Energy::writeBufferToH5()
{
	BinaryFile file(this->filename);
	file.addDatasetDouble(
		"time",
		this->times);
	file.addDatasetDouble(
		"energy",
		this->energies);
}

void Energy::writeBufferToDat()
{
	std::ofstream file;
	file.open(this->filename, std::ofstream::out);
	file << "Time\tEnergy\n";
	for (unsigned int i=0; i<this->energies.size(); i++) {
		file << this->times[i] << "\t" << this->energies[i] << "\n";
	}
	file.close();
}
