/* Particle.cpp
 * author: Christoph Froehner
 */

#include "Particle.h"

Particle::Particle()
{
	position 	= new double[3];
	force		= new double[3];
}

Particle::~Particle()
{
	delete position;
	delete force;
}

void Particle::move(double deviation[3])
{
	position[0] += deviation[0];
	position[1] += deviation[1];
	position[2] += deviation[2];
}

void Particle::addForce(double forceTerm[3])
{
	force[0] += forceTerm[0];
	force[1] += forceTerm[1];
	force[2] += forceTerm[2];
}
