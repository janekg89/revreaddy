/* Observable.h
 * 
 * Implement different observables objects, which will be 
 * recorded during simulation. 
 */

#ifndef __OBSERVABLE_H_INCLUDED__
#define __OBSERVABLE_H_INCLUDED__
#include <vector>
#include <string>
#include "Particle.h"

class Observable
{
	public:
		// number of timesteps between two recordings
		unsigned int recPeriod; 
		// number of timesteps between two writeBufferToFile
		unsigned int clearPeriod; 
		std::string filename;
		std::vector<unsigned long long> observedParticleIds;
		// buffer is declared by children
		virtual void record(
			std::vector<Particle>& activeParticles,
			double t); // write data to buffer
		virtual void writeBufferToFile(); // clear buffer
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
