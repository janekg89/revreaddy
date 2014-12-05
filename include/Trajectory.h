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
		// the struct particleTuple holds the type of a particle and its coordinates
		struct particleTuple
		{
			std::string particleType;
			std::array<double,3> particleCoordinates;
			unsigned long int particleTime;
		};
		// trajectory [time] [particles] [tuple[type][coordinates]]
		std::vector< std::vector< particleTuple > > trajectory;
		void record(std::vector<Particle> activeParticles, unsigned long int t);
		void writeBufferToFile();
	
		Trajectory();
		~Trajectory();
};

#endif
