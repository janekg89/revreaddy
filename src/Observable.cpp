/* Observable.cpp */

#include "Observable.h"

Observable::Observable() {}

Observable::~Observable() {}

void Observable::configure(World * world, Config * config) {}
void Observable::record(World * world, double t) {}

void Observable::writeBufferToFile() {}

// return the particle index within particles of particle with id id
int Observable::findParticleIndex(
	std::vector<Particle>& particles,
	unsigned long long id) 
{ 
	unsigned int max = particles.size() - 1;
	unsigned int min = 0;
	unsigned int mid = 0;
	while (max >= min) {
		mid = min + (max - min) / 2;
		if (particles[mid].uniqueId == id) {return mid;}
		else if (particles[mid].uniqueId < id) {min = mid + 1;}
		else {max = mid - 1;}
	}
	// std::cout << "id was not found\n"; //this has to be catched
	return -1;
}