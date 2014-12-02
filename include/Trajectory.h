/* Trajectory.h
 * 
 * Child class of Observable.
 * Records all particles positions.
 */

#ifndef __TRAJECTORY_H_INCLUDED__
#define __TRAJECTORY_H_INCLUDED__
#include "Observable.h"
#include "Particle.h"
#include <vector>
#include <array>
#include <fstream>
#include <string>

class Trajectory : public Observable
{
	public:
		// trajectory [time] [particles] [coordinate]
		std::vector< std::vector< std::array<double,3> > > trajectory;
		void record(std::vector<Particle> activeParticles, unsigned long int t);
		void writeBufferToFile();
	
		Trajectory();
		~Trajectory();
};

#endif
