/* MeanSquaredDisplacement.cpp */

#include "MeanSquaredDisplacement.h"

MeanSquaredDisplacement::MeanSquaredDisplacement(unsigned int id)
{
	this->recPeriod = 1;
	this->clearPeriod = 0;
	this->particleTypeId = id;
}

MeanSquaredDisplacement::~MeanSquaredDisplacement()
{
	this->recPeriod = 0;
}

void MeanSquaredDisplacement::record(
	std::vector<Particle>& activeParticles,
	double t)
{
	if (this->meanSquaredDisplacements.size() == 0) {
		for (int i=0; i<activeParticles.size(); i++) {
			if (activeParticles[i].typeId == this->particleTypeId) {
				positionTuple pt;
				pt.position = activeParticles[i].position;
				pt.boxCoordinates = activeParticles[i].boxCoordinates;
				this->startPoints.push_back(pt);
				this->startTime = t;
			}
		}
	}
	std::vector<double> displacement = {0.,0.,0.};
	std::vector<long>   relativeBoxCoordinates = {0,0,0};
	std::vector<double> squaredDisplacements;
	double squaredDisplacement = 0.;
	double mean = 0.;
	unsigned int index = 0;
	for (int i=0; i<activeParticles.size(); i++) {
		if (activeParticles[i].typeId == this->particleTypeId) {
			relativeBoxCoordinates[0] = activeParticles[i].boxCoordinates[0] 
			                          - startPoints[index].boxCoordinates[0];
			relativeBoxCoordinates[1] = activeParticles[i].boxCoordinates[1]
			                          - startPoints[index].boxCoordinates[1];
			relativeBoxCoordinates[2] = activeParticles[i].boxCoordinates[2]
			                          - startPoints[index].boxCoordinates[2];
			displacement[0] = (double) relativeBoxCoordinates[0] 
			                * this->boxsize
			                + activeParticles[i].position[0]
			                - startPoints[index].position[0];
			displacement[1] = (double) relativeBoxCoordinates[1]
			                * this->boxsize
			                + activeParticles[i].position[1]
			                - startPoints[index].position[1];
			displacement[2] = (double) relativeBoxCoordinates[2]
			                * this->boxsize
			                + activeParticles[i].position[2]
			                - startPoints[index].position[2];
			squaredDisplacement = displacement[0]*displacement[0]
			                    + displacement[1]*displacement[1]
			                    + displacement[2]*displacement[2];
			squaredDisplacements.push_back(squaredDisplacement);
			mean += squaredDisplacement;
			index += 1;
		}
	}
	mean /= (double) this->startPoints.size();
	this->meanSquaredDisplacements.push_back(mean);
	double standardDeviation = 0.;
	for (int j=0; j<squaredDisplacements.size(); j++) {
		standardDeviation += (squaredDisplacements[j] - mean)
		                   * (squaredDisplacements[j] - mean);
	}
	standardDeviation /= (double) (this->startPoints.size() - 1);
	double standardError = standardDeviation 
	                     / (double) this->startPoints.size();
	standardDeviation = sqrt(standardDeviation);
	standardError     = sqrt(standardError);
	this->standardDeviations.push_back(standardDeviation);
	this->standardErrors.push_back(standardError);
	this->time.push_back(t - this->startTime);
}

void MeanSquaredDisplacement::writeBufferToFile()
{
	std::ofstream file;
	file.open(this->filename);
	file << "Number of particles: " << startPoints.size() << "\n";
	file << "Time\tMean Squared Displacement\tStandardDeviation\tStandardError\n";
	for (int i=0; i<meanSquaredDisplacements.size(); i++) {
		file << time[i] << "\t";
		file << meanSquaredDisplacements[i] << "\t";
		file << standardDeviations[i] << "\t";
		file << standardErrors[i] << "\n";
	}
	file.close();
}
