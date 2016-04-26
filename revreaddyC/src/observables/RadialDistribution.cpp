/* RadialDistribution.cpp */

#include "RadialDistribution.h"

// receive vector<array<unsigned int>> representing a list of tuples
// these tuples have pairs of particleTypeIds which should be considered
// in rdf calculation.
RadialDistribution::RadialDistribution(unsigned long inRecPeriod, unsigned long inClearPeriod, std::vector<double>& range, std::vector< std::array<unsigned, 2> > considered, std::string inFilename) {
	this->recPeriod    = inRecPeriod;
	this->clearPeriod  = inClearPeriod;
	this->clearedAutomatically = false;
	isSetup = false;
	this->filename     = inFilename;
	observableTypeName = "RadialDistribution";
	this->numberOfBins = range.size() - 1;
	this->radialDistribution = gsl_histogram_alloc(this->numberOfBins);
	this->rangeOfBins  = range;
	const double * cRange = &range[0];
	gsl_histogram_set_ranges(this->radialDistribution, cRange, range.size());
	// calculate centers of bins
	double center;
	for (unsigned i=0; i < (rangeOfBins.size() - 1) ; i++) {
		center = 0.5 * ( rangeOfBins[i] + rangeOfBins[i+1] );
		this->binCenters.push_back(center);
		this->bins.push_back(0.);
	}
	this->consideredPairs = considered;
	// apply sorted convention
	for (auto&& pair : consideredPairs) {
		std::sort(pair.begin(), pair.end());
	}
	this->utils = new Utils();
}

RadialDistribution::~RadialDistribution() {
	gsl_histogram_free(this->radialDistribution);
	delete this->utils;
}

void RadialDistribution::configure(World * world, Config * config) {
	this->isPeriodic = config->isPeriodic;
	this->boxsize = config->boxsize;
}

/* Record the radial distribution already normalized 
 * correctly for the current timestep. */
void RadialDistribution::record(World * world, double t) {
	double radius = 0.;
	for (unsigned long i=0; i<world->particles.size(); i++) {
		for (unsigned long j=0; j<world->particles.size(); j++) {
			if (this->isInConsidered(
				world->particles[i].typeId,
				world->particles[j].typeId)) {
				if (i != j) {
					this->utils->getMinDistanceSquared(
						radius,
						world->particles[i].position,
						world->particles[j].position,
						this->isPeriodic,
						this->boxsize);
					radius = sqrt(radius);
					gsl_histogram_increment(this->radialDistribution, radius);
				}
			}
		}
	}
	// copy the hist to 'bins' while scaling every value correctly
	for (unsigned i=0; i<bins.size(); i++) {
		bins[i] += gsl_histogram_get(this->radialDistribution, i) / (binCenters[i] * binCenters[i]);
	}
	gsl_histogram_reset(this->radialDistribution);
}

/* Finds out if tuple (a,b) is in consideredPairs, this only
 * depends on the order of a and b. So you can only look at the
 * RDF from particle a to b solely. */
bool RadialDistribution::isInConsidered(unsigned a, unsigned b) {
	for (unsigned int k=0; k<this->consideredPairs.size(); k++) {
		if (this->consideredPairs[k][0] == a) {
			if (this->consideredPairs[k][1] == b) {
				return true;
			}
		}
	}
	return false;
}

void RadialDistribution::writeToH5() {
	H5::H5File file(this->filename.c_str(), H5F_ACC_TRUNC);
	createExtendibleDataset(file, "binCenters", this->binCenters);
	createExtendibleDataset(file, "bins", this->bins);
}

void RadialDistribution::writeToDat() {
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