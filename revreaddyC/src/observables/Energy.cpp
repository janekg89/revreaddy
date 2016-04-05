/* Energy.cpp */

#include "Energy.h"

Energy::Energy(unsigned inRecPeriod, unsigned inClearPeriod, std::string inFilename) {
	this->recPeriod = inRecPeriod;
	this->clearPeriod = inClearPeriod;
	this->clearedAutomatically = false;
	this->filename = inFilename;
	observableTypeName = "Energy";
	isSetup = false;
}

void Energy::record(World * world, double t) {
	this->energies.push_back(world->energy);
	this->times.push_back(t);
}

void Energy::writeToH5() {
	BinaryFile file(this->filename);
	file.addDatasetDouble(
		"time",
		this->times);
	file.addDatasetDouble(
		"energy",
		this->energies);
}

void Energy::writeToDat() {
	std::ofstream file;
	file.open(this->filename, std::ofstream::out);
	file << "Time\tEnergy\n";
	for (unsigned int i=0; i<this->energies.size(); i++) {
		file << this->times[i] << "\t" << this->energies[i] << "\n";
	}
	file.close();
}