/* RadialDistribution.cpp */

#include "RadialDistribution.h"

// receive vector<vector<unsigned int>> representing a list of tuples
// these tuples have pairs of particleTypeIds which should be considered
// in rdf calculation.
RadialDistribution::RadialDistribution(
	std::vector<double>& range,
	bool isPeriodic,
	double boxsize,
	std::vector< std::vector<unsigned int> > considered,
	std::string inFilename)
{
	this->recPeriod    = 1;
	this->clearPeriod  = 0;
	this->filename     = inFilename;
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
	World * world,
	double t)
{
	double radius = 0.;
	for (int i=0; i<world->activeParticles.size(); i++) {
		for (int j=0; j<world->activeParticles.size(); j++) {
			if (this->isInConsidered(
				world->activeParticles[i].typeId,
				world->activeParticles[j].typeId)) {
				if (i != j) {
					getMinDistanceSquared(
						radius,
						world->activeParticles[i].position,
						world->activeParticles[j].position,
						this->isPeriodic,
						this->boxsize);
					radius = sqrt(radius);
					gsl_histogram_increment(this->radialDistribution, radius);
				}
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
	// first determine the file extension
	unsigned int lastDot = this->filename.find_last_of(".");
	std::string extension = this->filename.substr(lastDot);
	if ( (extension == ".h5") || (extension == ".hdf5") ) {
		this->writeBufferToH5();
	}
	else if ( (extension == ".dat") || (extension == ".txt") ) {
		this->writeBufferToDat();
	}
	else {
		this->writeBufferToDat();
	}
}

void RadialDistribution::writeBufferToH5()
{
	BinaryFile file(this->filename);
	file.addDatasetDouble(
		"binCenters",
		this->binCenters);
	file.addDatasetDouble(
		"bins",
		this->bins);
}

void RadialDistribution::writeBufferToDat()
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