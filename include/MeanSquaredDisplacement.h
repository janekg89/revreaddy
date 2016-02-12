/* MeanSquaredDisplacement.h
 *
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

class MeanSquaredDisplacement : public Observable
{
	public:
		
		unsigned particleTypeId;
		// [particle] [x,y,z]
		struct positionTuple {
			std::vector<double> position;
			std::vector<long>   boxCoordinates;
		};
		std::vector< positionTuple > startPoints;
		double startTime;
		// [timeIndex]
		std::vector<double> time;
		std::vector<double> meanSquaredDisplacements;
		std::vector<double> standardDeviations;
		std::vector<double> standardErrors;
		std::vector<unsigned int> numberOfParticles;
		double boxsize;

		/* configure() saves the uniqueIds and initial positions of 
		 * the considered particles and saves the boxsize */
		void configure(World * world, Config * config);
		void record(World * world, double t);
		void writeBufferToFile();
		void writeBufferToH5();
		void writeBufferToDat();
		
		MeanSquaredDisplacement(
			unsigned long inRecPeriod,
			unsigned long inClearPeriod,
			unsigned int pTypeId,
			std::string inFilename);
		~MeanSquaredDisplacement();
};

#endif // __MEANSQUAREDDISPLACEMENT_H_INCLUDED__