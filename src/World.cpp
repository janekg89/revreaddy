/* World.cpp */

#include "World.h"

World::World()
{
	this->cumulativeRuntime     = 0.;
	this->energy                = 0.;
	this->oldEnergy             = 0.;	
	this->acceptProbDynamics    = 1.;
	this->acceptProbReactions   = 1.;
	this->acceptionsDynamics    = 0;
	this->rejectionsDynamics    = 0;
	this->uniqueIdCounter       = 0;
}

void World::addParticle(
	std::vector<double> initPos,
	unsigned int particleTypeId)
{
/*	if (particleTypeId >= config->particleTypes.size() ) {
		std::cout << "Error: The given particle type does not exist!\n"
		          << "Particle is not created" << std::endl;
		return;
	}*/
	Particle particle;
	if ( initPos.size() == 3 ) { particle.position  = initPos; }
	else {
		std::cout << "Error: Particles' initial position has dimension mismatch!\n" 
		          << "Particle will be placed at {0,0,0}" << std::endl;	
		particle.position = {0., 0., 0.};
	}
	particle.typeId = particleTypeId;
	particle.uniqueId = this->uniqueIdCounter;
	this->uniqueIdCounter += 1;
	this->particles.push_back(particle);//push_back copies arg into vec
}

void World::removeParticle(unsigned long int index)
{
	this->particles.erase(this->particles.begin() + index);
}

unsigned long World::getNumberOfParticles() { return this->particles.size(); }

std::vector<double> World::getPosition(unsigned long index) { return this->particles[index].position; }

void World::setPosition(unsigned long int index, std::vector<double> newPos)
{
	if (newPos.size() == 3) {
		this->particles[index].position[0] = newPos[0];
		this->particles[index].position[1] = newPos[1];
		this->particles[index].position[2] = newPos[2];
	}
	else {
		std::cout << "Error: New position has dimension mismatch!\n"
		          << "Particle remains at its old position" << std::endl;
	}
}

unsigned World::getTypeId(unsigned long index) {
	return this->particles[index].typeId;
}

void World::setTypeId(unsigned long index, unsigned typeId) 
{
/*	if (typeId >= config->particleTypes.size() ) {
		std::cout << "Error: The given particle type does not exist!\n"
		          << "Particle is not created" << std::endl;
		return;
	}*/
	this->particles[index].typeId = typeId;
}

void World::deleteAllParticles()
{
	 this->particles.erase(
		this->particles.begin(),
		this->particles.begin() + this->particles.size()
	);
}