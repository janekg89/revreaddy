/* ParticleNumbers.h
 * Count the number of particles of a specific type */

#ifndef __PARTICLENUMBERS_H_INCLUDED__
#define __PARTICLENUMBERS_H_INCLUDED__
#include <vector>
#include <string>
#include <fstream>
#include "World.h"
#include "Config.h"
#include "Observable.h"

class ParticleNumbers : public Observable
{
public:
	ParticleNumbers(
		unsigned long inRecPeriod,
		unsigned long inClearPeriod,
		std::string inFilename, 
		unsigned inParticleTypeId);
	void configure(World * world, Config * config);
	void record(World * world, double t);
	void writeBufferToFile();
	void writeBufferToH5();
	void writeBufferToDat();
private:
	unsigned particleTypeId;
	std::vector<unsigned long> particleNumbers;
	std::vector<double> time;
};

#endif //__PARTICLENUMBERS_H_INCLUDED__