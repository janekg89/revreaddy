/* Particle.h
 * author: Christoph Froehner
 *
 * This data structure will hold information
 * about a single particle. It will be constructed
 * for every particle in the simulation. */

#ifndef __PARTICLE_H_INCLUDED__
#define __PARTICLE_H_INCLUDED__
#include <iostream>
#include <vector>
class Particle
{
	public:
		unsigned long long uniqueId;
		unsigned int typeId; // determines potentials, diffConst and radius
		std::vector<double> position; // current position
		std::vector<long>   boxCoordinates; // id of box where particle is
		//unsigned int count; // how many timesteps can still be skipped
		/* how many timesteps in total were skipped. Determines 
		 * the distribution from which new position is drawn */
		//unsigned int skip;
		std::vector<double> cumulativeForce;

		Particle();

		void move(std::vector<double>& deviation);
		void addForce(std::vector<double>& force);
		void resetForce();
};

#endif
