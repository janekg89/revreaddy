/* Observable.cpp
 * 
 */

#include "Observable.h"

Observable::Observable()
{

}

Observable::~Observable()
{

}

void Observable::record(std::vector<Particle>& activeParticles, double t)
{

}

void Observable::writeBufferToFile()
{

}

// return the particle index within activeParticles of particle with id id
int findParticleIndex(
	std::vector<Particle>& activeParticles,
	unsigned long long id) 
{ 
	unsigned int max = activeParticles.size() - 1;
	unsigned int min = 0;
	unsigned int mid = 0;
	while (max >= min) {
		mid = min + (max - min) / 2;
		if (activeParticles[mid].uniqueId == id) {return mid;}
		else if (activeParticles[mid].uniqueId < id) {min = mid + 1;}
		else {max = mid - 1;}
	}
	// std::cout << "id was not found\n"; //this has to be catched
	return -1;
} 
