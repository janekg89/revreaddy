/* MeanSquaredDisplacement.cpp */

#include "MeanSquaredDisplacement.h"

MeanSquaredDisplacement::MeanSquaredDisplacement(
	std::vector<Particle>& activeParticles,
	unsigned int pTypeId,
	double time,
	double boxsize,
	std::string inFilename)
{
	this->recPeriod = 1;
	this->clearPeriod = 0;
	this->filename = inFilename;
	this->particleTypeId = pTypeId;
	this->startTime = time;
	this->boxsize   = boxsize;
	/* add uniqueIds of all particles of type 
	 * "particleTypeId" to observedParticleIds */
	for (unsigned int i=0; i<activeParticles.size(); i++) {
		if (activeParticles[i].typeId == this->particleTypeId) {
			this->observedParticleIds.push_back( activeParticles[i].uniqueId );
			positionTuple pt;
			pt.position = activeParticles[i].position;
			pt.boxCoordinates = activeParticles[i].boxCoordinates;
			this->startPoints.push_back(pt);
		}
	}
}

MeanSquaredDisplacement::~MeanSquaredDisplacement()
{
}

void MeanSquaredDisplacement::record(
	World * world,
	double t)
{
	std::vector<double> displacement = {0.,0.,0.};
	std::vector<long>   relativeBoxCoordinates = {0,0,0};
	std::vector<double> squaredDisplacements;
	double squaredDisplacement = 0.;
	double mean = 0.;
	int index = 0;
	unsigned int stillExisting = 0;
	for (unsigned i=0; i < this->observedParticleIds.size(); i++) {
		// check if particle still exists. if not: index = -1
		index = this->findParticleIndex(
			world->activeParticles,
			observedParticleIds[i]);
		if (index != -1) {
			stillExisting += 1;
			relativeBoxCoordinates[0] =world->activeParticles[index].boxCoordinates[0]
			                          - startPoints[i].boxCoordinates[0];
			relativeBoxCoordinates[1] =world->activeParticles[index].boxCoordinates[1]
			                          - startPoints[i].boxCoordinates[1];
			relativeBoxCoordinates[2] =world->activeParticles[index].boxCoordinates[2]
			                          - startPoints[i].boxCoordinates[2];
			displacement[0] = (double) relativeBoxCoordinates[0] 
			                * this->boxsize
			                + world->activeParticles[index].position[0]
			                - startPoints[i].position[0];
			displacement[1] = (double) relativeBoxCoordinates[1]
			                * this->boxsize
			                + world->activeParticles[index].position[1]
			                - startPoints[i].position[1];
			displacement[2] = (double) relativeBoxCoordinates[2]
			                * this->boxsize
			                + world->activeParticles[index].position[2]
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
	for (unsigned j=0; j<squaredDisplacements.size(); j++) {
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
	this->time.push_back(t - this->startTime);
}

void MeanSquaredDisplacement::writeBufferToFile()
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

void MeanSquaredDisplacement::writeBufferToH5()
{
	BinaryFile file(this->filename);
	file.addDatasetDouble(
		"time",
		this->time);
	file.addDatasetDouble(
		"meanSquaredDisplacement",
		this->meanSquaredDisplacements);
	file.addDatasetDouble(
		"standardDeviation",
		this->standardDeviations);
	file.addDatasetDouble(
		"standardError",
		this->standardErrors);
	file.addDatasetUnsignedInt(
		"numberOfParticles",
		this->numberOfParticles);
}

void MeanSquaredDisplacement::writeBufferToDat()
{
	std::ofstream file;
	file.open(this->filename);
	file << "Time\tMeanSquaredDisplacement\tStandardDeviation\t"
			"StandardError\tNumberOfParticles\n";
	for (unsigned i=0; i<meanSquaredDisplacements.size(); i++) {
		file << time[i] << "\t";
		file << meanSquaredDisplacements[i] << "\t";
		file << standardDeviations[i] << "\t";
		file << standardErrors[i] << "\t";
		file << numberOfParticles[i] << "\n";
	}
	file.close();
}