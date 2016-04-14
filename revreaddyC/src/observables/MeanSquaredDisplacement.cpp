/* MeanSquaredDisplacement.cpp */

#include "MeanSquaredDisplacement.h"

MeanSquaredDisplacement::MeanSquaredDisplacement(unsigned long inRecPeriod, unsigned long inClearPeriod, unsigned pTypeId, std::string inFilename) {
	this->recPeriod = inRecPeriod;
	this->clearedAutomatically = false;
	this->clearPeriod = inClearPeriod;
	this->filename = inFilename;
	observableTypeName = "MeanSquaredDisplacement";
	this->particleTypeId = pTypeId;
	isSetup = false;
}

void MeanSquaredDisplacement::configure(World * world, Config * config) {
	/* keep track of system parameters that might have changed */
	this->boxsize = config->boxsize;
	this->startTime = world->cumulativeRuntime;
}

void MeanSquaredDisplacement::setup(World * world, Config * config) {
	/* add uniqueIds of all particles of type "particleTypeId" to observedParticleIds */
	for (auto i=0; i<world->particles.size(); ++i) {
		if (world->particles[i].typeId == this->particleTypeId) {
			this->observedParticleIds.push_back( world->particles[i].uniqueId );
			positionTuple pt;
			pt.position = world->particles[i].position;
			pt.boxCoordinates = world->particles[i].boxCoordinates;
			this->startPoints.push_back(pt);
		}
	}
}

void MeanSquaredDisplacement::record(World * world, double t) {
	std::vector<double> displacement = {0.,0.,0.};
	std::vector<long> relativeBoxCoordinates = {0,0,0};
	std::vector<double> squaredDisplacements;
	double squaredDisplacement = 0.;
	double mean = 0.;
	int index = 0;
	unsigned int stillExisting = 0;
	for (auto i=0; i < this->observedParticleIds.size(); i++) {
		// check if particle still exists. if not: index = -1
		index = this->findParticleIndex(
			world->particles,
			observedParticleIds[i]);
		if (index != -1) {
			stillExisting += 1;
			relativeBoxCoordinates[0] =world->particles[index].boxCoordinates[0]
			                          - startPoints[i].boxCoordinates[0];
			relativeBoxCoordinates[1] =world->particles[index].boxCoordinates[1]
			                          - startPoints[i].boxCoordinates[1];
			relativeBoxCoordinates[2] =world->particles[index].boxCoordinates[2]
			                          - startPoints[i].boxCoordinates[2];
			displacement[0] = (double) relativeBoxCoordinates[0] 
			                * this->boxsize
			                + world->particles[index].position[0]
			                - startPoints[i].position[0];
			displacement[1] = (double) relativeBoxCoordinates[1]
			                * this->boxsize
			                + world->particles[index].position[1]
			                - startPoints[i].position[1];
			displacement[2] = (double) relativeBoxCoordinates[2]
			                * this->boxsize
			                + world->particles[index].position[2]
			                - startPoints[i].position[2];
			squaredDisplacement = displacement[0]*displacement[0]
			                    + displacement[1]*displacement[1]
			                    + displacement[2]*displacement[2];
			squaredDisplacements.push_back(squaredDisplacement);
			mean += squaredDisplacement;
		}
	}
	// if no more particles still exist, stop here. nothing is recorded
	if (stillExisting == 0) {return;}
	this->numberOfParticles.push_back(stillExisting);
	mean /= (double) stillExisting;
	this->meanSquaredDisplacements.push_back(mean);
	double standardDeviation = 0.;
	for (auto j=0; j<squaredDisplacements.size(); j++) {
		standardDeviation += (squaredDisplacements[j] - mean)
		                   * (squaredDisplacements[j] - mean);
	}
	standardDeviation /= (double) (stillExisting - 1);
	double standardError = standardDeviation 
	                     / (double) stillExisting; 
	standardDeviation = sqrt(standardDeviation);
	standardError     = sqrt(standardError);
	this->standardDeviations.push_back(standardDeviation);
	this->standardErrors.push_back(standardError);
	this->times.push_back(t - this->startTime);
}

void MeanSquaredDisplacement::writeToH5() {
	H5::H5File file(this->filename.c_str(), H5F_ACC_TRUNC);
	createExtendibleDataset(file, "times", this->times);
	createExtendibleDataset(file, "meanSquaredDisplacements", this->meanSquaredDisplacements);
	createExtendibleDataset(file, "standardDeviations", this->standardDeviations);
	createExtendibleDataset(file, "standardErrors", this->standardErrors);
	createExtendibleDataset(file, "numberOfParticles", this->numberOfParticles);
	this->times.clear();
	this->meanSquaredDisplacements.clear();
	this->standardDeviations.clear();
	this->standardErrors.clear();
	this->numberOfParticles.clear();
}

void MeanSquaredDisplacement::writeToDat() {
	std::ofstream file;
	file.open(this->filename);
	file << "times\tmeanSquaredDisplacement\tstandardDeviation\t"
			"standardError\tnumberOfParticles\n";
	for (unsigned i=0; i<meanSquaredDisplacements.size(); i++) {
		file << times[i] << "\t";
		file << meanSquaredDisplacements[i] << "\t";
		file << standardDeviations[i] << "\t";
		file << standardErrors[i] << "\t";
		file << numberOfParticles[i] << "\n";
	}
	file.close();
	this->times.clear();
	this->meanSquaredDisplacements.clear();
	this->standardDeviations.clear();
	this->standardErrors.clear();
	this->numberOfParticles.clear();
}