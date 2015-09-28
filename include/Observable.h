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
		// buffer is declared by children
		virtual void record(
			World * world,
			double t); // write data to buffer
		// TODO in the future only use flushing to avoid
		// bugs occuring when writing data, that has already
		// been written.
		/* write to file but keep the buffer in memory */
		virtual void writeBufferToFile();
		/* write to file and clear the buffer in memory */
		//virtual void flushBufferToFile();
		int findParticleIndex(
			std::vector<Particle>& activeParticles,
			unsigned long long id);
		
		/* define a 'virtual' variable to access Trajectory's trajectory 
		 * from Simulationn, which only knows about Observable. */
		std::vector< std::vector<double> > trajectory;
		
		Observable();
		virtual ~Observable();		
};

#endif // __OBSERVABLE_H_INCLUDED__
