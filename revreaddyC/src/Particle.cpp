#include "Particle.h"
Particle::Particle() {
	this->uniqueId = 0;
	this->position = {0.,0.,0.};
	this->position.shrink_to_fit();
	this->boxCoordinates = {0,0,0};
	this->boxCoordinates.shrink_to_fit();
	this->cumulativeForce = {0.,0.,0.};
	this->cumulativeForce.shrink_to_fit();
	this->typeId = 0;
}

void Particle::move(std::vector<double>& deviation) {
	this->position[0] += deviation[0];
	this->position[1] += deviation[1];
	this->position[2] += deviation[2];
}

void Particle::addForce(std::vector<double>& force) {
	this->cumulativeForce[0] += force[0];
	this->cumulativeForce[1] += force[1];
	this->cumulativeForce[2] += force[2];
}

void Particle::resetForce() {
	this->cumulativeForce = {0., 0., 0.};
}