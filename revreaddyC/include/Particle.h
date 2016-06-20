/* Particle.h
 * This data structure will hold information
 * about a single particle. It will be constructed
 * for every particle in the simulation. */

#ifndef __PARTICLE_H_INCLUDED__
#define __PARTICLE_H_INCLUDED__
#include <iostream>
#include <vector>
class Particle {
public:
	unsigned long long uniqueId;
	unsigned int typeId; // determines potentials, diffConst and radius
	std::vector<double> position; // current position
	std::vector<long> boxCoordinates; // id of box where particle is
	std::vector<double> cumulativeForce;

	Particle();

	void move(std::vector<double>& deviation);
	void addForce(std::vector<double>& force);
	void resetForce();
};

#endif //__PARTICLE_H_INCLUDED__