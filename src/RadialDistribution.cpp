/* RadialDistribution.cpp */

#include "RadialDistribution.h"

// TODO receive vector<vector<unsigned int>> representing a list of tuples
// these tuples have pairs of particleTypeIds which should be considered
// in rdf calculation.
RadialDistribution::RadialDistribution(
	std::vector<double>& range,
	bool isPeriodic,
	double boxsize,
	std::vector< std::vector<unsigned int> > considered)
{
	this->recPeriod    = 1;
	this->clearPeriod  = 0;
	this->isPeriodic   = isPeriodic;
	this->boxsize      = boxsize;
	this->numberOfBins = range.size() - 1;
	this->radialDistribution = gsl_histogram_alloc(this->numberOfBins);
	this->rangeOfBins  = range;
	const double * cRange = &range[0];
	gsl_histogram_set_ranges(this->radialDistribution, cRange, range.size());
	// calculate centers of bins
	double center;
	for (int i=0; i < (rangeOfBins.size() - 1) ; i++) {
		center = 0.5 * ( rangeOfBins[i] + rangeOfBins[i+1] );
		this->binCenters.push_back(center);
		this->bins.push_back(0.);
	}
	this->consideredPairs = considered;
}

RadialDistribution::~RadialDistribution()
{
	gsl_histogram_free(this->radialDistribution);
}

/* Record the radial distribution already normalized 
 * correctly for the current timestep.
 */
void RadialDistribution::record(
	std::vector<Particle>& activeParticles,
	double t)
{
	double radius = 0.;
	for (int i=0; i<activeParticles.size(); i++) {
		for (int j=i+1; j<activeParticles.size(); j++) {
			if (this->isInConsidered(
				activeParticles[i].typeId,
				activeParticles[j].typeId)) {
					getMinDistanceSquared(
						radius,
						activeParticles[i].position,
						activeParticles[j].position,
						this->isPeriodic,
						this->boxsize);
					radius = sqrt(radius);
					gsl_histogram_increment(this->radialDistribution, radius);
			}
		}
	}
	// copy the hist to 'bins' while scaling every value correctly
	for (int i=0; i<bins.size(); i++) {
		bins[i] += gsl_histogram_get(this->radialDistribution, i) 
		           / (binCenters[i] * binCenters[i]);
	}
	gsl_histogram_reset(this->radialDistribution);
}

/* Finds out if tuple (a,b) is in consideredPairs*/
bool RadialDistribution::isInConsidered(unsigned int a, unsigned int b)
{
	if (a <= b) {
		for (unsigned int k=0; k<this->consideredPairs.size(); k++) {
			if (this->consideredPairs[k][0] == a) {
				if (this->consideredPairs[k][1] == b) {
					return true;
				}
			}
		}
		return false;
	}
	else {
		for (unsigned int k=0; k<this->consideredPairs.size(); k++) {
			if (this->consideredPairs[k][0] == b) {
				if (this->consideredPairs[k][1] == a) {
					return true;
				}
			}
		}
		return false;
	}
}

void RadialDistribution::writeBufferToFile()
{
	std::ofstream file;
	file.open(this->filename, std::ofstream::out);
	for (auto&& center : this->binCenters) {
		file << center << "\t";
	}
	file << "\n";
	for (auto&& bin : this->bins) {
		file << bin  << "\t";
	}
	file.close();
}
