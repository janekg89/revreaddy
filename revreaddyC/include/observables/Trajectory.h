/* Trajectory.h - Child class of Observable. Records all particles' positions */
#ifndef __TRAJECTORY_H_INCLUDED__
#define __TRAJECTORY_H_INCLUDED__
#include "World.h"
#include "Config.h"
#include "Particle.h"
#include <vector>
#include <fstream>
#include <string>
#include "Observable.h"

class Trajectory : public Observable {
public:
	void record(World * world, double t);
	/* Write trajectory to file in append-like fashion. 
	 * Trajectory buffer is cleared afterwards. */
	void writeToDat(); // xyz format
	void writeToH5(); // so far no h5 implementation. not efficient.

	Trajectory(unsigned long inRecPeriod, unsigned long inClearPeriod, std::string inFilename);

private:
	/* A snapshot contains the system particles' information at a single point in time. */
	struct Snapshot {
		std::vector<double> positionsX;
		std::vector<double> positionsY;
		std::vector<double> positionsZ;
		std::vector<double> forcesX;
		std::vector<double> forcesY;
		std::vector<double> forcesZ;
		std::vector<unsigned int> typeIds;
		std::vector<unsigned long long> uniqueIds;
		double time;
		unsigned long timeIncrement;
		unsigned long numberOfParticles;
	};
	/* The data buffer where all information is stored. Is flushed upon writeToFile(). */
	std::vector< Snapshot > trajectory;
	/* tells which snapshot has to be recorded to buffer next */
	unsigned long nextRecordedTimeIncrement;
};

#endif // __TRAJECTORY_H_INCLUDED__