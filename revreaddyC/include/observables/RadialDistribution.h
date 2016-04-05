/* RadialDistribution.h
 * This class is a child of Observable.
 * Calculate the Radial Distribution Function cumulative over many timesteps.
 * Make use of the GSL libraries to calculate the RDF as a histogram. */

#ifndef __RADIALDISTRIBUTION_H_INCLUDED__
#define __RADIALDISTRIBUTION_H_INCLUDED__
#include <vector>
#include <array>
#include <fstream>
#include <string>
#include <gsl/gsl_histogram.h>
#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
#include "World.h"
#include "Config.h"
#include "Particle.h"
#include "utils.h"
#include "Observable.h"

class RadialDistribution : public Observable {
public:
	void configure(World * world, Config * config);
	void record(World * world, double t);
	/* Check if (a,b) is in consideredPairs.*/
	bool isInConsidered(unsigned a, unsigned b);
	void writeToH5();
	void writeToDat();

	RadialDistribution(
		unsigned long inRecPeriod,
		unsigned long inClearPeriod,
		std::vector<double>& range,
		std::vector< std::array<unsigned, 2> > considered,
		std::string inFilename);
	~RadialDistribution();

private:
	gsl_histogram * radialDistribution;
	size_t numberOfBins;
	std::vector<double> rangeOfBins;
	std::vector<double> binCenters;
	/* bins is basically a copy of the histogram to easy rescale each
	 * bin individually. */
	std::vector<double> bins;
	/* consideredPairs is a list of tuples of particleTypeIds
	 * CONVENTION: the tuples MUST be ordered.
	 * Correct: (0,1) , (3,5), (2,2) 
	 * Not correct: (1,0), (4,1)
	 * This is because isInConsidered() searches for tuples (a,b)
	 * where a<=b.*/
	std::vector< std::array<unsigned, 2> > consideredPairs;
	bool isPeriodic;
	double boxsize;
	Utils * utils;
};

#endif // __RADIALDISTRIBUTION_H_INCLUDED__