/* Fusion.cpp
 * @todo this has to be reworked. especially reactionDistance shall be used.
 * simplify configuration */
#include "Fusion.h"
#include <math.h>
#include <armadillo>


Fusion::Fusion(
	std::string inName,
	std::vector<unsigned> inForwardTypes,
	std::vector<unsigned> inBackwardTypes,
	double inForwardRate,
	double inBackwardRate,
	double inReactionDistance)
{
	this->name = inName;
	this->forwardTypes = inForwardTypes;
	this->backwardTypes = inBackwardTypes;
	this->forwardRate = inForwardRate;
	this->backwardRate = inBackwardRate;
	this->reactionDistance = inReactionDistance;
	this->type = "Fusion";
	this->utils = new Utils();
}

Fusion::~Fusion()
{
	delete this->utils;
}

void Fusion::configure(
	std::vector< std::shared_ptr<ParticleInteraction> > inInteractions,
	//double inInversePartition,
	double inMaxDistr,
	double inMeanDistr,
	double inInverseTemperature,
	double inRadiusA,
	double inRadiusB,
	bool inIsPeriodic,
	double inBoxsize,
	double inWeightA,
	double inWeightB)
{
	this->interactions = inInteractions;
	//this->inversePartition = inInversePartition;
	this->maxDistr = inMaxDistr;
	this->meanDistr = inMeanDistr;
	this->inverseTemperature = inInverseTemperature;
	this->radiusA = inRadiusA;
	this->radiusB = inRadiusB;
	//this->weightA = pow(radiusA, 3.) / (pow(radiusA, 3.)+pow(radiusB, 3.));
	//this->weightB = pow(radiusB, 3.) / (pow(radiusA, 3.)+pow(radiusB, 3.));
	weightA = inWeightA;
	weightB = inWeightB;
	this->isPeriodic = inIsPeriodic;
	this->boxsize = inBoxsize;
	radiiSum = radiusA + radiusB;
}

double Fusion::performForward(
	std::vector<unsigned long> particleIndices,
	double timestep,
	World * world,
	Random * random)
{
	double forwardProb = this->forwardRate * timestep;
	double u = random->uniform();
	if ( u < forwardProb ) {
		/* reaction occurs */
		unsigned long indexI = particleIndices[0];
		unsigned long indexJ = particleIndices[1];
		/* r_ij is the vector from particle i to j */
		std::vector<double> r_ij = {0.,0.,0.};
		this->utils->getMinDistanceVector(
			r_ij,
			world->particles[indexI].position,
			world->particles[indexJ].position,
			this->isPeriodic,
			this->boxsize);
		double weight;
		if (world->particles[indexI].typeId == forwardTypes[0]) {
			weight = this->weightB;
		}
		else {
			weight = this->weightA;
		}
		/* position is the point in between particles i and j,
		 * weighted according to their masses
		 * pos = r_A + (r_B - r_A) * weightB, with r_ij = r_B - r_A */
		std::vector<double> position = {0.,0.,0.};
		position[0] = world->particles[indexI].position[0] + weight * r_ij[0];
		position[1] = world->particles[indexI].position[1] + weight * r_ij[1]; 
		position[2] = world->particles[indexI].position[2] + weight * r_ij[2];

		/* SUPER IMPORTANT: delete the particle with the
         * higher index first or the indexing in particles
         * gets mixed up and you delete wrong particles */
		if (world->useFractional) {
			auto alpha = world->alpha;
			if ( indexI > indexJ ) {
				world->removeParticleAndIncrements(indexI);
				world->removeParticleAndIncrements(indexJ);
			}
			else {
				world->removeParticleAndIncrements(indexJ);
				world->removeParticleAndIncrements(indexI);
			}
			world->addParticleAndIncrements(position, backwardTypes[0], random, maxTime, timestep, diffC, alpha);
		} else {
			if ( indexI > indexJ ) {
				world->removeParticle(indexI);
				world->removeParticle(indexJ);
			}
			else {
				world->removeParticle(indexJ);
				world->removeParticle(indexI);
			}
			world->addParticle(position, this->backwardTypes[0]);
		}
		double distance = sqrt(r_ij[0]*r_ij[0]+r_ij[1]*r_ij[1]+r_ij[2]*r_ij[2]);
        /* increment reaction counter */
        world->forwardReactionCounter[name]++;
		/* f(d) is part of the forward proposal probability */
		return distribution( distance );
	}
	else {/* nothing happens */ return 1.;}
}

double Fusion::performBackward(
	std::vector<unsigned long> particleIndices,
	double timestep,
	World * world,
	Random * random)
{
	double backwardProb = this->backwardRate * timestep;
	double u = random->uniform();
	if ( u < backwardProb ) {
		/* reaction occurs */
		unsigned long index = particleIndices[0];
		std::vector<double> positionC = world->particles[index].position;
		/* distance = [0, reactionRadiiSum] */
		double distance = this->randomFromDistribution(random);
		/* orientation is a random unit vector, i.e. with length = 1 */
		std::vector<double> orientation = random->normal3D();
		double norm = sqrt(
			  orientation[0]*orientation[0]
			+ orientation[1]*orientation[1]
			+ orientation[2]*orientation[2]);
		orientation[0] /= norm;
		orientation[1] /= norm;
		orientation[2] /= norm;
		/* determine new positions for particles A and B,
		 * the distance from positionC weighted with the mass
		 * of the other particle */
		std::vector<double> positionA = positionC;
		std::vector<double> positionB = positionC;
		positionA[0] += this->weightB * distance * orientation[0];
		positionA[1] += this->weightB * distance * orientation[1];
		positionA[2] += this->weightB * distance * orientation[2];
		positionB[0] -= this->weightA * distance * orientation[0];
		positionB[1] -= this->weightA * distance * orientation[1];
		positionB[2] -= this->weightA * distance * orientation[2];
		/* remove C, add A and B */
		if (world->useFractional) {
			auto alpha = world->alpha;
			world->removeParticleAndIncrements(index);
			world->addParticleAndIncrements(positionA, forwardTypes[0], random, maxTime, timestep, diffA, alpha);
			world->addParticleAndIncrements(positionB, forwardTypes[1], random, maxTime, timestep, diffB, alpha);
		} else {
			world->removeParticle(index);
			world->addParticle(positionA, this->forwardTypes[0]);
			world->addParticle(positionB, this->forwardTypes[1]);
		}
        /* increment reaction counter */
        world->backwardReactionCounter[name]++;
		/* 1./f(d) is part of the backward proposal probability */
		return 1./distribution(distance);
	}
	else {/* nothing happens */ return 1.;}
}

double Fusion::distribution(double x) {
	double result = 0.;
	for (unsigned i=0; i<this->interactions.size(); i++) {
		result += this->interactions[i]->calculateEnergy(
			x * x,// * this->reactionRadiiSum * this->reactionRadiiSum, // distance (x * R_AB) squared
			this->radiiSum * this->radiiSum); // sum of radii squared
	}
	result = exp(- this->inverseTemperature * result);
	//result *= this->inversePartition;
	return result;
}

double Fusion::randomFromDistribution(Random * random) {
	unsigned it = 0;
	double x;
	double y;
	while ( it < 100 ) {
		x = random->uniform() * this->reactionDistance;
		y = this->maxDistr * random->uniform();
		if ( y < this->distribution(x) ) { return x; }
		else { it += 1; }
	}
	LOG_INFO("uniformFromDistribution left with mean value.")
	return this->meanDistr;
}

void Fusion::configureFractional(unsigned long inMaxTime, double inDiffA, double inDiffB, double inDiffC) {
	this->maxTime = inMaxTime;
	this->diffA = inDiffA;
	this->diffB = inDiffB;
	this->diffC = inDiffC;
}
