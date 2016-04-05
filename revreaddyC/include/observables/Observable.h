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

class Observable {
public:
	// decide if the observable should be written during runtime or manually
	bool clearedAutomatically;
	// some observables require setup that may only happen once in contrast to configure
	bool isSetup;
	// number of timesteps between two recordings
	unsigned long int recPeriod; 
	// number of recordings between two writeBufferToFile
	// not yet used
	unsigned long int clearPeriod; 
	std::string filename;
	std::string observableTypeName;
	std::vector<unsigned long long> observedParticleIds;
	/* configure sets up the observables' parameters that could have
	 * changed between construction of the observable and the start
	 * of the simulation */
	virtual void configure(World * world, Config * config);
	/* setup cares about the initialization direclty befire the run. In
	 * contrast to configure, this must only happen once, e.g. saving
	 * uniqueIds of considered particles. */
	virtual void setup(World * world, Config * config);
	/* record data from World and write to the buffer */
	virtual void record(World * world, double t);
	/* write to file and clear the buffer */
	void writeToFile();
	virtual void writeToH5();
	virtual void writeToDat();
	long int findParticleIndex(std::vector<Particle>& particles, unsigned long long id);

	Observable();
	virtual ~Observable();		
};

#endif // __OBSERVABLE_H_INCLUDED__
