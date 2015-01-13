/* Observable.h
 * 
 * Implement different observables objects, which will be 
 * recorded during simulation. 
 */

#ifndef __OBSERVABLE_H_INCLUDED__
#define __OBSERVABLE_H_INCLUDED__
#include <vector>
#include "Particle.h"

class Observable
{
	public:
		int recPeriod; // number of timesteps between two recordings
		int clearPeriod; // number of timesteps between two writeBufferToFile
		// buffer is declared by children
		virtual void record(
			std::vector<Particle> activeParticles,
			unsigned long int t); // write data to buffer
		virtual void writeBufferToFile(); // clear buffer
		
		Observable();
		virtual ~Observable();		
};
#endif // __OBSERVABLE_H_INCLUDED__
