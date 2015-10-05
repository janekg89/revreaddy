/* Fusion.cpp */
#include "Fusion.h"
#include <math.h>

Fusion::Fusion(
	std::string inName,
	std::vector<unsigned> inForwardTypes,
	std::vector<unsigned> inBackwardTypes,
	double inForwardRate,
	double inBackwardRate,
	Random * inRandom)
{
	this->name = inName;
	this->forwardTypes = inForwardTypes;
	this->backwardTypes = inBackwardTypes;
	this->forwardRate = inForwardRate;
	this->backwardRate = inBackwardRate;
	this->random = inRandom;
	this->type = "Fusion";
	this->isConfigured = false;
}

void Fusion::configure(
	std::vector< std::shared_ptr<ParticleInteraction> > inInteractions,
	double inInversePartition,
	double inMaxDistr,
	double inRadiiSum,
	double inMeanDistr,
	double inInverseTemperature)
{
	this->interactions = inInteractions;
	this->inversePartition = inInversePartition;
	this->maxDistr = inMaxDistr;
	this->radiiSum = inRadiiSum;
	this->meanDistr = inMeanDistr;
	this->inverseTemperature = inInverseTemperature;
	this->isConfigured = true;
}

double Fusion::performForward(
	std::vector<unsigned long> particleIndices,
	World * world,
	double timestep)
{
	double forwardProb = this->forwardRate * timestep;
	double u = this->random->uniform();
	if ( u < forwardProb ) {
		/* reaction occurs */
		unsigned long indexI = particleIndices[0];
		unsigned long indexJ = particleIndices[1];
		/* position is the point in between particles i and j */
		std::vector<double> position = {0.,0.,0.};
		position[0] = 0.5 * (
			  world->activeParticles[indexI].position[0]
			+ world->activeParticles[indexJ].position[0]);
		position[1] = 0.5 * (
			  world->activeParticles[indexI].position[1]
			+ world->activeParticles[indexJ].position[1]);
		position[2] = 0.5 * (
			  world->activeParticles[indexI].position[2]
			+ world->activeParticles[indexJ].position[2]);
		/* SUPER IMPORTANT: delete the particle with the 
		 * higher index first or the indexing in activeParticles
		 * gets mixed up and you delete wrong particles */
		if ( indexI > indexJ ) {
			world->removeParticle(indexI);
			world->removeParticle(indexJ);
		}
		else {
			world->removeParticle(indexJ);
			world->removeParticle(indexI);
		}
		world->addParticle(
			position,
			this->backwardTypes[0]);
		return 1.;
	}
	else {/* nothing happens */ return 1.;}
}

double Fusion::performBackward(
	std::vector<unsigned long> particleIndices,
	World * world,
	double timestep)
{
	double backwardProb = this->backwardRate * timestep;
	double u = this->random->uniform();
	if ( u < backwardProb ) {
		/* reaction occurs */
		unsigned long index = particleIndices[0];
		std::vector<double> positionC = world->activeParticles[index].position;
		double distance = this->radiiSum * this->uniformFromDistribution();
		// TODO orientation is non-uniform w.r.t. phi and theta.
		// make it uniform!
		/* orientation is a random unit vector, i.e. with length = 1 */
		std::vector<double> orientation = {0.,0.,0.};
		orientation[0] = this->random->uniform() - 0.5;
		orientation[1] = this->random->uniform() - 0.5; 
		orientation[2] = this->random->uniform() - 0.5; 
		double norm = sqrt(
			  orientation[0]*orientation[0]
			+ orientation[1]*orientation[1]
			+ orientation[2]*orientation[2]);
		orientation[0] /= norm;
		orientation[1] /= norm;
		orientation[2] /= norm;
		/* determine new positions for particles A and B */
		std::vector<double> positionA = positionC;
		std::vector<double> positionB = positionC;
		positionA[0] += 0.5 * distance * orientation[0];
		positionA[1] += 0.5 * distance * orientation[1];
		positionA[2] += 0.5 * distance * orientation[2];
		positionB[0] -= 0.5 * distance * orientation[0];
		positionB[1] -= 0.5 * distance * orientation[1];
		positionB[2] -= 0.5 * distance * orientation[2];
		/* remove C, add A and B */
		world->removeParticle(index);
		world->addParticle(positionA, this->forwardTypes[0]);
		world->addParticle(positionB, this->forwardTypes[1]);
		return 1.;
	}
	else {/* nothing happens */ return 1.;}
}

double Fusion::distribution(double x)
{
	double result = 0.;
	for (unsigned i=0; this->interactions.size(); i++) {
		result += this->interactions[i]->calculateEnergy(
			x * x * this->radiiSum * this->radiiSum, // distance (x * R_AB) squared
			this->radiiSum * this->radiiSum); // sum of radii squared
	}
	result = exp(- this->inverseTemperature * result);
	result *= x * this->inversePartition;
	return result;
}

double Fusion::uniformFromDistribution()
{
	unsigned it = 0;
	double x = 0.;
	double y = 0.;
	while ( it < 100 ) {
		x = this->random->uniform();
		y = this->maxDistr * this->random->uniform();
		if ( y < this->distribution(x) ) { return x; }
		else { it += 1; }
	}
	return this->meanDistr;
}