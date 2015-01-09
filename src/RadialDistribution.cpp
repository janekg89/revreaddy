/* RadialDistribution.cpp */

#include "RadialDistribution.h"

RadialDistribution::RadialDistribution(std::vector<double> range)
{
	this->numberOfBins = range.size() - 1;
	this->radialDistribution = gsl_histogram_alloc(this->numberOfBins);
	this->rangeOfBins = range;
	const double * cRange = &range[0];
	gsl_histogram_set_ranges(this->radialDistribution, cRange, range.size());
	// calculate centers of bins
	double center;
	for (int i=0; i < (rangeOfBins.size() - 1) ; i++) {
		center = 0.5 * ( rangeOfBins[i] + rangeOfBins[i+1] );
		this->binCenters.push_back(center);
		this->bins.push_back(0.);
	}
}

RadialDistribution::~RadialDistribution()
{
	gsl_histogram_free(this->radialDistribution);
}


void RadialDistribution::record(std::vector<Particle> activeParticles, unsigned long int t)
/* Record the radial distribution already normalized 
 * correctly for the current timestep.
 */
{
	double radius;
	std::array<double,3> r_ij;
	for (int i=0; i<activeParticles.size(); i++) {
		for (int j=i+1; j<activeParticles.size(); j++) {
			r_ij = getMinDistanceVector(
				activeParticles[i].position,
				activeParticles[j].position,
				this->isPeriodic,
				this->boxsize 
			);
			radius = r_ij[0]*r_ij[0] + r_ij[1]*r_ij[1] + r_ij[2]*r_ij[2];
			radius = sqrt(radius);
			gsl_histogram_increment(this->radialDistribution, radius);
		}
	}
	// copy the hist to 'bins' while scaling every value correctly
	for (int i=0; i<bins.size(); i++) {
		bins[i] += gsl_histogram_get(this->radialDistribution, i) / (binCenters[i] * binCenters[i]);
	}
	gsl_histogram_reset(this->radialDistribution);
}

void RadialDistribution::writeBufferToFile()
{
	std::ofstream file;
	file.open("radialdistribution.dat");
	for (auto&& center : this->binCenters) {
		file << center << "\t";
	}
	file << "\n";
	for (auto&& bin : this->bins) {
		file << bin  << "\t";
	}
	file.close();
}
