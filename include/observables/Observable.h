/* Observable.h
 * 
 * Implement different observables objects, which will be 
 * recorded during simulation. 
 */

#ifndef __OBSERVABLE_H_INCLUDED__
#define __OBSERVABLE_H_INCLUDED__
#include <vector>
#include <string>
#include <fstream>
#include "World.h"
#include "Config.h"
#include "Particle.h"
#include "BinaryFile.h"

class Observable
{
public:
	// number of timesteps between two recordings
	unsigned long int recPeriod; 
	// number of recordings between two writeBufferToFile
	// not yet used
	unsigned long int clearPeriod; 
	std::string filename;
	std::vector<unsigned long long> observedParticleIds;
	/* configure sets up the observables' parameters that could have
	 * changed between construction of the observable and the start
	 * of the simulation */
	virtual void configure(World * world, Config * config);
	virtual void record(World * world, double t); // write data to buffer
	// TODO in the future only use flushing to avoid
	// bugs occuring when writing data, that has already
	// been written.
	/* write to file but keep the buffer in memory */
	virtual void writeBufferToFile();
	/* write to file and clear the buffer in memory */
	//virtual void flushBufferToFile();
	int findParticleIndex(
		std::vector<Particle>& particles,
		unsigned long long id);

	Observable();
	virtual ~Observable();		
};

#endif // __OBSERVABLE_H_INCLUDED__
