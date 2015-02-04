/* MeanSquaredDisplacement.h
 *
 * Child class of Observable. Records the MeanSquaredDisplacement of a
 * specific particle type */

#ifndef __MEANSQUAREDDISPLACEMENT_H_INCLUDED__
#define __MEANSQUAREDDISPLACEMENT_H_INCLUDED__
#include "Observable.h"
#include "Particle.h"
#include <vector>
#include <fstream>
#include <math.h>

class MeanSquaredDisplacement : public Observable
{
	public:
		
		unsigned int particleTypeId;
		// [particle] [x,y,z]
		struct positionTuple {
			std::vector<double> position;
			std::vector<long>   boxCoordinates;
		};
		std::vector< positionTuple > startPoints;
		// [timeIndex]
		std::vector<double> meanSquaredDisplacements;
		std::vector<double> standardDeviations;
		std::vector<double> standardErrors;
		double boxsize;

		void record(
			std::vector<Particle>& activeParticles,
			unsigned long int t);
		void writeBufferToFile();
		
		MeanSquaredDisplacement(unsigned int id);
		~MeanSquaredDisplacement();
};

#endif // __MEANSQUAREDDISPLACEMENT_H_INCLUDED__
