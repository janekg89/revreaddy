/* ProbabilityDensity.cpp */

#include "ProbabilityDensity.h"

ProbabilityDensity::ProbabilityDensity(
	unsigned long inRecPeriod,
	unsigned long inClearPeriod,
	std::string inFilename,
	unsigned particleTypeId,
	std::vector<double>& range,
	unsigned coord)
{
	this->recPeriod = inRecPeriod;
	this->clearPeriod = inClearPeriod;
	this->filename = inFilename;
	if (coord >= 3) {
		throw Exception("The given coordinate can only be 0, 1 or 2");
	}
	this->coordinate = coord;
	this->rangeOfBins = range;
	this->numberOfBins = range.size() - 1;
	this-> probabilityDensity = gsl_histogram_alloc(this->numberOfBins);
	const double * cRange = &range[0];
	gsl_histogram_set_ranges(this->probabilityDensity, cRange, range.size());
	/* calculate centers of bins */
	double center;
	for (unsigned i=0; i < (rangeOfBins.size() - 1) ; i++) {
		center = 0.5 * (rangeOfBins[i] + rangeOfBins[i+1] );
		this->binCenters.push_back(center);
		this->bins.push_back(0.);
	}
	this->particleTypeId = particleTypeId;
}

ProbabilityDensity::~ProbabilityDensity()
{
}

void ProbabilityDensity::configure(World * world, Config * config)
{
	/* add uniqueIds of all particles of type 
	 * "particleTypeId" to observedParticleIds */
	for (unsigned i=0; i<world->particles.size(); i++) {
		if (world->particles[i].typeId == this->particleTypeId) {
			this->observedParticleIds.push_back( world->particles[i].uniqueId );
		}
	}
}

void ProbabilityDensity::record(World * world, double t)
{
	int index = 0;
	for (unsigned i=0; i < this->observedParticleIds.size(); i++) {
		// check if particle still exists. if not: index = -1
		index = this->findParticleIndex(
			world->particles,
			observedParticleIds[i]);
		if (index != -1) {
			gsl_histogram_increment(
				this->probabilityDensity,
				world->particles[index].position[this->coordinate]);
		}
	}
	// copy the hist to "bins" and reset the gsl hist
	for (unsigned i=0; i<bins.size(); i++) {
		bins[i] += gsl_histogram_get(this->probabilityDensity, i);
	}
	gsl_histogram_reset(this->probabilityDensity);
}

void ProbabilityDensity::writeBufferToFile()
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

void ProbabilityDensity::writeBufferToH5()
{
	BinaryFile file(this->filename);
	file.addDatasetDouble(
		"binCenters",
		this->binCenters);
	file.addDatasetDouble(
		"bins",
		this->bins);
}

void ProbabilityDensity::writeBufferToDat()
{
	std::ofstream file;
	file.open(this->filename, std::ofstream::out);
	for (auto&& center : this->binCenters) {
		file << center << "\t";
	}
	file << "\n";
	for (auto&& bin : this->bins) {
		file << bin << "\t";
	}
	file.close();
}
