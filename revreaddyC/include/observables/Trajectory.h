/* Trajectory.h
 * Child class of Observable. Records all particles' positions */

#ifndef __TRAJECTORY_H_INCLUDED__
#define __TRAJECTORY_H_INCLUDED__
#include "World.h"
#include "Config.h"
#include "Particle.h"
#include <vector>
#include <fstream>
#include <string>
#include "Observable.h"

class Trajectory : public Observable 
{
public:
	void configure(World * world, Config * config);
	void record(World * world, double t);
	void writeBufferToFile();
	void writeBufferToH5();
	void writeBufferToDat();

	Trajectory(unsigned long inRecPeriod, unsigned long inClearPeriod, std::string inFilename);
	~Trajectory();

private:
	/* the struct particleTuple holds the type of a particle 
	 * and its coordinates */
	struct particleTuple {
		unsigned int particleTypeId;
		std::vector<double> particleCoordinates;
		double particleTime;
	};
	/* trajectory has the following shape:
	 * [relativeTime] [particles] [tuple[type][coordinates][absTime]] */
	std::vector< std::vector< particleTuple > > trajectory;
};

#endif // __TRAJECTORY_H_INCLUDED__