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
		/* consideredPairs is a list of tuples of particleTypeIds
		 * CONVENTION: the tuples MUST be ordered.
		 * Correct: (0,1) , (3,5), (2,2) 
		 * Not correct: (1,0), (4,1)
		 * This is because isInConsidered() searches for tuples (a,b)
		 * where a<=b.*/
		std::vector< std::vector<unsigned int> > consideredPairs;
		bool isPeriodic;
		double boxsize;

		void record(
			std::vector<Particle>& activeParticles,
			double t);
		/* Check if (a,b) is in consideredPairs.*/
		bool isInConsidered(unsigned int a, unsigned int b);
		void writeBufferToFile();

		RadialDistribution(
			std::vector<double>& range,
			bool isPeriodic,
			double boxsize,
			std::vector< std::vector<unsigned int> > considered,
			std::string inFilename);
		~RadialDistribution();
};
#endif // __RADIALDISTRIBUTION_H_INCLUDED__
