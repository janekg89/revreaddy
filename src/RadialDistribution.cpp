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
{
	double radius;
	std::array<double,3> r_ij;
	for (int i=0; i<activeParticles.size(); i++) {
		for (int j=i; j<activeParticles.size(); j++) {
			r_ij = this->sim->getMinDistance( activeParticles[i].position, activeParticles[j].position );
			radius = r_ij[0]*r_ij[0] + r_ij[1]*r_ij[1] + r_ij[2]*r_ij[2];
			radius = sqrt(radius);
			gsl_histogram_increment(this->radialDistribution, radius);
		}
	}
}

void RadialDistribution::writeBufferToFile()
{
	std::ofstream file;
	file.open("radialdistribution.dat");
	for (auto&& val : this->rangeOfBins) {
		file << val << "\t";
	}
	file << "\n";
	for (int i=0; i<numberOfBins; i++) {
		file << gsl_histogram_get(this->radialDistribution, i) << "\t";
	}
	file << "0";
}
