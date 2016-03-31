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
	/* Write trajectory to file in append-like fashion. Trajectory buffer is cleared
	 * afterwards. */
	void writeBufferToH5();
	void writeBufferToDat(); // xyz format

	Trajectory(unsigned long inRecPeriod, unsigned long inClearPeriod, std::string inFilename);
	~Trajectory();

private:
	/* the struct particleTuple holds the type of a particle 
	 * and its coordinates */
	struct particleTuple {
		unsigned int particleTypeId;
		unsigned long long particleUniqueId;
		std::vector<double> particlePosition;
		std::vector<double> particleForce;
		double particleTime;
	};
	/* trajectory has the following shape:
	 * [relativeTime] [particles] [tuple[type][coordinates][absTime]] */
	std::vector< std::vector< particleTuple > > trajectory;
	/* tells which group (time index) has to be added next */
	unsigned long nextWrittenTimeIncrement;
};

#endif // __TRAJECTORY_H_INCLUDED__