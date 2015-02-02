/* Particle.h
 * author: Christoph Froehner
 *
 * This data structure will hold information
 * about a single particle. It will be constructed
 * for every particle in the simulation.
 * 
 * Particles' typeId governs which potentials are used
 * and which reactions are performed. The convention is:
 * 0 - None (will behave like a soft particle, but should
 *     not be used though)
 * 1 - Lennard-Jones particle
 * 2 - soft particle
 * ... all coming numbers are reserved for different user
 *     specific soft particles.
 */

#ifndef __PARTICLE_H_INCLUDED__
#define __PARTICLE_H_INCLUDED__
#include <iostream>
#include <vector>
class Particle
{
	public:
		unsigned int typeId; // determines potentials
		std::vector<double> position; // current position
		std::vector<long>   boxCoordinates; // id of box where particle is
		std::vector<double> oldPosition; // position when forces were calc'd
		std::vector<long>   oldBoxCoordinates;
		unsigned int count; // how many timesteps can still be skipped
		/* how many timesteps in total were skipped. Determines 
		 * the distribution from which new position is drawn */
		unsigned int skip;
		double radius;	// size of the particle
		/* an individual diffusionConstant, typically determined by the type */
		double diffusionConstant; 
		/* the cumulative force for the current timestep. */
		std::vector<double> cumulativeForce;
		std::vector<double> oldForce;

		Particle();
		~Particle();
		void move(std::vector<double> deviation);
		void addForce(std::vector<double> force);
		void resetForce();
};

#endif
