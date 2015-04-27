/* ProbabilityDensity.h 
 * 
 * Record a histogram of positions in x or y or z. Can observe one or several
 * particles. Vanished particles are not considered anymore. */

#ifndef __PROBABILITYDENSITY_H_INCLUDED__
#define __PROBABILITYDENSITY_H_INCLUDED__
#include "Observable.h"
#include "Particle.h"
#include <vector>
#include <fstream>
#include <gsl/gsl_histogram.h>

class ProbabilityDensity : public Observable
{
	public:
		gsl_histogram * probabilityDensity;
		/* The coordinate (x,y,z) = (0,1,2), of which to record
		 * the probability density. */
		unsigned int coordinate;
		size_t numberOfBins;
		std::vector<double> rangeOfBins;
		std::vector<double> binCenters;
		/* bins is basically a copy of the histogram to easy rescale each
		 * bin individually. */
		std::vector<double> bins;
		unsigned int particleTypeId;

		void record(
			std::vector<Particle>& activeParticles,
			double t);
		void writeBufferToFile();
		void writeBufferToH5();
		void writeBufferToDat();

		ProbabilityDensity(
			std::vector<Particle>& activeParticles,
			unsigned int pTypeId,
			std::vector<double>& range,
			unsigned int coord,
			std::string inFilename);
		~ProbabilityDensity();
};

#endif // __PROBABILITYDENSITY_H_INCLUDED__
