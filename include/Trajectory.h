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
#include <fstream>
#include <string>

class Trajectory : public Observable
{
	public:
		/* the struct particleTuple holds the type of a particle 
		 * and its coordinates */
		struct particleTuple
		{
			unsigned int particleTypeId;
			std::vector<double> particleCoordinates;
			double particleTime;
		};
		/* trajectory has the following shape:
		 * [relativeTime] [particles] [tuple[type][coordinates][absTime]] */
		std::vector< std::vector< particleTuple > > trajectory;

		void record(
			World * world,
			double t);
		void writeBufferToFile();
	
		Trajectory(std::string inFilename);
		~Trajectory();
};

#endif // __TRAJECTORY_H_INCLUDED__