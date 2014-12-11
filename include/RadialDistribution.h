/* RadialDistribution.h
 * 
 * This class is a child of Observable.
 * Calculate the Radial Distribution Function cumulative over many timesteps.
 * Make use of the GSL libraries to calculate the RDF as a histogram.
 *
 * Histogram is written to file in the following format
 *  binEdge0	binEdge1	binEdge2 ... binEdgeN
 *  binContent0	binContent1	binContent2 ... binContentN-1
 */

#ifndef __RADIALDISTRIBUTION_H_INCLUDED__
#define __RADIALDISTRIBUTION_H_INCLUDED__
#include "Observable.h"
#include "Particle.h"
#include "Simulation.h" // use getMinDistance(r1, r2)
#include <array>
#include <vector>
#include <fstream>
#include <string>
#include <gsl/gsl_histogram.h>
#define _USE_MATH_DEFINES
#include <cmath>

class RadialDistribution : public Observable
{
	public:
		gsl_histogram * radialDistribution;
		size_t numberOfBins;
		std::vector<double> rangeOfBins;
		std::vector<double> binCenters;
		/* bins is basically a copy of the histogram to easy rescale each
		 * bin individually. */
		std::vector<double> bins;
		std::string particleType; // the type for which to calculate the RDF

		Simulation * sim;//pointer to the Simulation that owns this observable

		void record(std::vector<Particle> activeParticles, unsigned long int t);
		void writeBufferToFile();
		void setRange(std::vector<double> range);

		RadialDistribution(size_t bins, Simulation * simulation);
		~RadialDistribution();
};
#endif // __RADIALDISTRIBUTION_H_INCLUDED__
