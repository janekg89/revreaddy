/* Particle.cpp
 * author: Christoph Froehner
 */

#include "Particle.h"

Particle::Particle()
{
	this->position 	= new double[3];
	this->force		= new double[3];
	this->count		= 0;
	this->diffusionConstant	= 1.;
	this->numberOfTimestep	= 0;
	this->radius	= 1.;
	this->skip		= 0;
}

Particle::~Particle()
{
	delete this->position;
	delete this->force;
}

void Particle::move(double deviation[3])
{
	this->position[0] += deviation[0];
	this->position[1] += deviation[1];
	this->position[2] += deviation[2];
}

void Particle::addForce(double forceTerm[3])
{
	this->force[0] += forceTerm[0];
	this->force[1] += forceTerm[1];
	this->force[2] += forceTerm[2];
}

void Particle::resetForce()
{
	this->force[0] = 0.;
	this->force[1] = 0.;
	this->force[2] = 0.;
}
