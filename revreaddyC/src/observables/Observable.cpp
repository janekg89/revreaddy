/* Observable.cpp */

#include "Observable.h"

Observable::Observable() {}

Observable::~Observable() {}

void Observable::configure(World * world, Config * config) {}
void Observable::record(World * world, double t) {}

void Observable::writeToFile() {
	// first determine the file extension
	unsigned int lastDot = this->filename.find_last_of(".");
	std::string extension = this->filename.substr(lastDot);
	if ( (extension == ".h5") || (extension == ".hdf5") ) {
		writeToH5();
	}
	else if ( (extension == ".dat") || (extension == ".txt") ) {
		writeToDat();
	}
	else {
		// write to hdf5 format per default
		writeToH5();
	}
}

void Observable::writeToH5() {
	LOG_TRACE("Enter Observable::writeToH5")
}

void Observable::writeToDat() {
	LOG_TRACE("Enter Observable::writeToDat")
}

// return the particle index within particles of particle with id id
long int Observable::findParticleIndex(std::vector<Particle>& particles, unsigned long long id) { 
	signed long max = particles.size() - 1;
	signed long min = 0;
	signed long mid = 0;
	while (max >= min) {
		mid = min + (max - min) / 2;
		if (particles[mid].uniqueId == id) {return mid;}
		else if (particles[mid].uniqueId < id) {min = mid + 1;}
		else {max = mid - 1;}
	}
	// std::cout << "id was not found\n"; //this has to be catched
	return -1;
}