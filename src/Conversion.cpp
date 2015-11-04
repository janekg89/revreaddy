/* Conversion.cpp */

#include "Conversion.h"

//#define __DEBUG__
#ifndef __DEBUG__
#define print(a)  
#endif
#ifdef __DEBUG__
#define print(a) std::cout << a << std::endl;
#undef __DEBUG__
#endif

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
	this->type = "Conversion";
}

Conversion::~Conversion()
{

}

double Conversion::performForward(
	std::vector<unsigned long int> particleIndices,
	World * world,
	double timestep)
{
	print("Enter Conversion performForward")
	/* approximation to Poisson probability */
	double forwardProb = this->forwardRate * timestep;
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
		return 1.;
	}
	else {/* nothing happens */ return 1.;}
}

double Conversion::performBackward(
	std::vector<unsigned long int> particleIndices,
	World * world,
	double timestep)
{
	print("Enter Conversion performBackward")
	/* approximation to Poisson probability */
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
		return 1.;
	}
	else {/* nothing happens */ return 1.;}
}