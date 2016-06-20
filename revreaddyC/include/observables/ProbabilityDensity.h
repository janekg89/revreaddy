/* ProbabilityDensity.h 
 * Record a histogram of positions in x or y or z. Can observe one or several
 * particles. Vanished particles are not considered anymore. */

#ifndef __PROBABILITYDENSITY_H_INCLUDED__
#define __PROBABILITYDENSITY_H_INCLUDED__
#include <vector>
#include <fstream>
#include <gsl/gsl_histogram.h>
#include "Observable.h"
#include "World.h"
#include "Config.h"
#include "Particle.h"
#include "Exception.h"

class ProbabilityDensity : public Observable {
public:
	/* setup saves the uniqueIds f considered particles. */
	void setup(World * world, Config * config);
	void record(World * world, double t);
	void writeToH5();
	void writeToDat();

	ProbabilityDensity(unsigned long inRecPeriod, unsigned long inClearPeriod, std::string inFilename, unsigned particleTypeId,	std::vector<double>& range,	unsigned coord);
private:
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
};

#endif // __PROBABILITYDENSITY_H_INCLUDED__
