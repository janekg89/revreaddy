/* MeanSquaredDisplacement.h
 * Child class of Observable. Records the MeanSquaredDisplacement of a
 * specific particle type */

#ifndef __MEANSQUAREDDISPLACEMENT_H_INCLUDED__
#define __MEANSQUAREDDISPLACEMENT_H_INCLUDED__
#include <vector>
#include <fstream>
#include <math.h>
#include "World.h"
#include "Config.h"
#include "Particle.h"
#include "Observable.h"

class MeanSquaredDisplacement : public Observable {
public:
	/* setup() saves the uniqueIds and initial positions of 
	 * the considered particles and configure() saves the boxsize */
	void setup(World * world, Config * config);
	void configure(World * world, Config * config);
	void record(World * world, double t);
	void writeToH5();
	void writeToDat();
	
	MeanSquaredDisplacement(unsigned long inRecPeriod, unsigned long inClearPeriod, unsigned int pTypeId, std::string inFilename);

private:
	unsigned particleTypeId;
	// [particle] [x,y,z]
	struct positionTuple {
		std::vector<double> position;
		std::vector<long>   boxCoordinates;
	};
	std::vector< positionTuple > startPoints;
	double startTime;
	std::vector<double> times;
	std::vector<double> meanSquaredDisplacements;
	std::vector<double> standardDeviations;
	std::vector<double> standardErrors;
	std::vector<unsigned int> numberOfParticles;
	double boxsize;
};

#endif // __MEANSQUAREDDISPLACEMENT_H_INCLUDED__