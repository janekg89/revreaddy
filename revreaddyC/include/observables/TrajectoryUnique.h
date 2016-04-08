/* TrajectoryUnique.h - Record positions of a set of predefined particles via their uniqueId.
 * I.e. if they vanish, they are not considered anymore (position is written as [0,0,0]).
 * More important, the list of particles never grows, which enables efficient storage. 
 * Currently consider all particles that were in the system initially. */
#ifndef __TRAJECTORYUNIQUE_H_INCLUDED__
#define __TRAJECTORYUNIQUE_H_INCLUDED__
#include <vector>
#include <map>
#include <string>
#include <boost/multi_array.hpp>
#include "Observable.h" 
#include "World.h"
#include "Config.h"
#include "Particle.h"

class TrajectoryUnique : public Observable {
public:
	TrajectoryUnique(unsigned long inRecPeriod, unsigned long inClearPeriod, std::string inFilename);

	void setup(World * world, Config * config);
	void record(World * world, double t);
	void writeToDat();
	void writeToH5();
private:
	unsigned long numberOfParticles;
	std::vector<unsigned long long> uniqueIds;
	std::vector<unsigned> typeIds;
	/* system information at a single point in time. */
	struct Snapshot {
		boost::multi_array<double,2> positions; // N x 3
		boost::multi_array<double,2> forces; // N x 3
		std::vector<bool> exists; // N, only important if particles vanish
		double time;
	};
	/* The data buffer */
	std::vector< Snapshot > trajectory;

	/* This is to check during record() if the particle with uniqueIds[i]
	 * does still exist. This value then is stillExists[i]. */
	std::map<unsigned long, bool> stillExists;

	// This is called when the given file already exists
	// ATTENTION: this will even append to existing trajectories that might have
	// been created with another system setup. This will lead to undetermined behavior.
	void appendToH5();
	// This is called when the given file does not exist
	void writeToNewH5();
	// copy data chunk from trajectory to 'real' arrays
	void bufferTimeDependentData(boost::multi_array<double,3>& positions, boost::multi_array<double,3>& forces, boost::multi_array<bool,2>& exists, std::vector<double>& times);
};

#endif //__TRAJECTORYUNIQUE_H_INCLUDED__