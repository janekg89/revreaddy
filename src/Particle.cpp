/* Particle.cpp
 * author: Christoph Froehner
 */

#include "Particle.h"
Particle::Particle()
{
	this->position 			= {0.,0.,0.};
	this->cumulativeForce 	= {0.,0.,0.};
	this->position.shrink_to_fit();
	this->cumulativeForce.shrink_to_fit();
	this->oldPosition       = {0.,0.,0.};
	this->oldPosition.shrink_to_fit();
	this->type				= "soft";
	this->diffusionConstant	= 1.;
	this->radius            = 1.;
	this->skip              = 0;
	this->count	            = 0;
	this->singleEnergy      = 0.;
}

Particle::~Particle()
{
}

void Particle::move(std::vector<double> deviation)
{
	this->position[0] += deviation[0];
	this->position[1] += deviation[1];
	this->position[2] += deviation[2];
}

void Particle::addForce(std::vector<double> force)
{
	this->cumulativeForce[0] += force[0];
	this->cumulativeForce[1] += force[1];
	this->cumulativeForce[2] += force[2];
}

void Particle::resetForce()
{
	this->cumulativeForce = {0., 0., 0.};
}
