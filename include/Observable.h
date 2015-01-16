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
		int recPeriod; // number of timesteps between two recordings
		int clearPeriod; // number of timesteps between two writeBufferToFile
		std::string filename;
		// buffer is declared by children
		virtual void record(
			std::vector<Particle> activeParticles,
			unsigned long int t); // write data to buffer
		virtual void writeBufferToFile(); // clear buffer
		
		/* define a 'virtual' variable to access Trajectory's trajectory 
		 * from Simulationn, which only knows about Observable. */
		std::vector< std::vector<double> > trajectory;
		
		Observable();
		virtual ~Observable();		
};
#endif // __OBSERVABLE_H_INCLUDED__
