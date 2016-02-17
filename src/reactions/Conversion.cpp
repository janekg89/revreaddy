/* Conversion.cpp */

#include "Conversion.h"

Conversion::Conversion(
	std::string inName,
	std::vector<unsigned int> inForwardTypes,
	std::vector<unsigned int> inBackwardTypes,
	double inForwardRate,
	double inBackwardRate)
{
	this->name = inName;
	this->forwardTypes = inForwardTypes;
	this->backwardTypes = inBackwardTypes;
	this->forwardRate = inForwardRate;
	this->backwardRate = inBackwardRate;
	this->reactionDistance = 0.;
	this->type = "Conversion";
}

Conversion::~Conversion() {}

double Conversion::performForward(
	std::vector<unsigned long int> particleIndices,
	double timestep,
	World * world,
	Random * random)
{
	/* approximation to Poisson probability */
	double forwardProb = this->forwardRate * timestep;
	double u = random->uniform();
	if ( u < forwardProb ) {
		/* reaction occurs */
		unsigned long index = particleIndices[0];
		std::vector<double> position = {0.,0.,0.};
		position[0] = world->particles[index].position[0];
		position[1] = world->particles[index].position[1];
		position[2] = world->particles[index].position[2];
		world->removeParticle(index);
		world->addParticle(
			position,
			this->backwardTypes[0]);
		return 1.;
	}
	else {/* nothing happens */ return 1.;}
}

double Conversion::performBackward(
	std::vector<unsigned long int> particleIndices,
	double timestep,
	World * world,
	Random * random)
{
	/* approximation to Poisson probability */
	double backwardProb = this->backwardRate * timestep;
	double u = random->uniform();
	if ( u < backwardProb ) {
		/* reaction occurs */
		unsigned long index = particleIndices[0];
		std::vector<double> position = {0.,0.,0.};
		position[0] = world->particles[index].position[0];
		position[1] = world->particles[index].position[1];
		position[2] = world->particles[index].position[2];
		world->removeParticle(index);
		world->addParticle(
			position,
			this->forwardTypes[0]);
		return 1.;
	}
	else {/* nothing happens */ return 1.;}
}