/* Conversion.cpp */

#include "Conversion.h"

Conversion::Conversion(
	std::string inName,
	std::vector<unsigned int> inForwardTypes,
	std::vector<unsigned int> inBackwardTypes,
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
}

Conversion::~Conversion()
{

}

double Conversion::performForward(
	std::vector<unsigned long int> particleIndices,
	World * world,
	double timestep)
{
	/* approximation to Poisson probability */
	double forwardProb = this->forwardRate * timestep;
	double backwardProb = this->backwardRate * timestep;
	double u = this->random->uniform();
	if ( u < forwardProb ) {
		/* reaction occurs */
		unsigned long index = particleIndices[0];
		std::vector<double> position = {0.,0.,0.};
		position[0] = world->activeParticles[index].position[0];
		position[1] = world->activeParticles[index].position[1];
		position[2] = world->activeParticles[index].position[2];
		world->removeParticle(index);
		world->addParticle(
			position,
			this->backwardTypes[0]);
		return (backwardProb / forwardProb);
	}
	else {/* nothing happens*/ return 1.;}
}

double Conversion::performBackward(
	std::vector<unsigned long int> particleIndices,
	World * world,
	double timestep)
{
	/* approximation to Poisson probability */
	double forwardProb = this->forwardRate * timestep;
	double backwardProb = this->backwardRate * timestep;
	double u = this->random->uniform();
	if ( u < backwardProb ) {
		/* reaction occurs */
		unsigned long index = particleIndices[0];
		std::vector<double> position = {0.,0.,0.};
		position[0] = world->activeParticles[index].position[0];
		position[1] = world->activeParticles[index].position[1];
		position[2] = world->activeParticles[index].position[2];
		world->removeParticle(index);
		world->addParticle(
			position,
			this->forwardTypes[0]);
		return (forwardProb / backwardProb);
	}
	else {/* nothing happens*/ return 1.;}
}