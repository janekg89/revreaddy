/* ProbabilityDensity.cpp */

#include "ProbabilityDensity.h"

ProbabilityDensity::ProbabilityDensity(unsigned long inRecPeriod, unsigned long inClearPeriod, std::string inFilename, unsigned particleTypeId, std::vector<double>& range, unsigned coord) {
	this->recPeriod = inRecPeriod;
	this->clearedAutomatically = false;
	this->clearPeriod = inClearPeriod;
	this->filename = inFilename;
	observableTypeName = "ProbabilityDensity";
	isSetup = false;
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

void ProbabilityDensity::setup(World * world, Config * config) {
	/* add uniqueIds of all particles of type 
	 * "particleTypeId" to observedParticleIds */
	for (auto i=0; i<world->particles.size(); i++) {
		if (world->particles[i].typeId == this->particleTypeId) {
			this->observedParticleIds.push_back( world->particles[i].uniqueId );
		}
	}
}

void ProbabilityDensity::record(World * world, double t) {
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

void ProbabilityDensity::writeToH5() {
	H5::H5File file(this->filename.c_str(), H5F_ACC_TRUNC);
	createExtendibleDataset(file, "binCenters", this->binCenters);
	createExtendibleDataset(file, "bins", this->bins);
}

void ProbabilityDensity::writeToDat() {
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
