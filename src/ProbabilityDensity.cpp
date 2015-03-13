/* ProbabilityDensity.cpp */

#include "ProbabilityDensity.h"

ProbabilityDensity::ProbabilityDensity(
	std::vector<Particle>& activeParticles,
	unsigned int pTypeId,
	std::vector<double>& range,
	unsigned int coord)
{
	if (coord < 3) {this->coordinate = coord;}
	else {coord = 0;}
	this->rangeOfBins = range;
	this->numberOfBins = range.size() - 1;
	this-> probabilityDensity = gsl_histogram_alloc(this->numberOfBins);
	const double * cRange = &range[0];
	gsl_histogram_set_ranges(this->probabilityDensity, cRange, range.size());
	/* calculate centers of bins */
	double center;
	for (int i=0; i < (rangeOfBins.size() - 1) ; i++) {
		center = 0.5 * (rangeOfBins[i] + rangeOfBins[i+1] );
		this->binCenters.push_back(center);
		this->bins.push_back(0.);
	}
	this->particleTypeId = pTypeId;
	/* add uniqueIds of all particles of type 
	 * "particleTypeId" to observedParticleIds */
	for (unsigned int i=0; i<activeParticles.size(); i++) {
		if (activeParticles[i].typeId == this->particleTypeId) {
			this->observedParticleIds.push_back( activeParticles[i].uniqueId );
		}
	}
}

ProbabilityDensity::~ProbabilityDensity()
{
}

void ProbabilityDensity::record(
	std::vector<Particle>& activeParticles,
	double t)
{
	int index = 0;
	for (int i=0; i < this->observedParticleIds.size(); i++) {
		// check if particle still exists. if not: index = -1
		index = this->findParticleIndex(
			activeParticles,
			observedParticleIds[i]);
		if (index != -1) {
			gsl_histogram_increment(
				this->probabilityDensity,
				activeParticles[index].position[this->coordinate]);
		}
	}
	// copy the hist to "bins" and reset the gsl hist
	for (int i=0; i<bins.size(); i++) {
		bins[i] += gsl_histogram_get(this->probabilityDensity, i);
	}
	gsl_histogram_reset(this->probabilityDensity);
}

void ProbabilityDensity::writeBufferToFile()
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
