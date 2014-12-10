/* RadialDistribution.cpp */

#include "RadialDistribution.h"

RadialDistribution::RadialDistribution(size_t bins, Simulation * simulation)
{
	this->radialDistribution = gsl_histogram_alloc(bins);
	this->numberOfBins = bins;
	this->sim = simulation;
}

RadialDistribution::~RadialDistribution()
{
	gsl_histogram_free(this->radialDistribution);
}

void RadialDistribution::setRange(std::vector<double> range)
{
	this->rangeOfBins = range;
	const double * cRange = &range[0];
	gsl_histogram_set_ranges(this->radialDistribution, cRange, range.size());
}

void RadialDistribution::record(std::vector<Particle> activeParticles, unsigned long int t)
/* Record the radial distribution already weighted 
 * correctly for the current timestep.
 */
{
	double radius;
	double weight;
	std::array<double,3> r_ij;
	for (int i=0; i<activeParticles.size(); i++) {
		for (int j=i+1; j<activeParticles.size(); j++) {
			r_ij = this->sim->getMinDistance( activeParticles[i].position, activeParticles[j].position );
			radius = r_ij[0]*r_ij[0] + r_ij[1]*r_ij[1] + r_ij[2]*r_ij[2];
			weight = this->sim->boxsize * this->sim->boxsize * this->sim->boxsize;
			weight /= 4. * M_PI * radius * (double) activeParticles.size();
			radius = sqrt(radius);
			gsl_histogram_accumulate(this->radialDistribution, radius, weight);
		}
	}
}

void RadialDistribution::writeBufferToFile()
{
	double center;
	std::ofstream file;
	file.open("radialdistribution.dat");
	// calculate centers of bins and write them to file
	for (int i=0; i < (rangeOfBins.size() - 1) ; i++) {
		center = 0.5 * ( rangeOfBins[i] + rangeOfBins[i+1] );
		file << center << "\t";
	}
	file << "\n";
	for (int i=0; i<numberOfBins; i++) {
		file << gsl_histogram_get(this->radialDistribution, i) << "\t";
	}
	file.close();
}
