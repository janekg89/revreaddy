/* Particle.h
 * author: Christoph Froehner
 *
 * This data structure will hold information
 * about a single particle. It will be constructed
 * for every particle in the simulation.
 *
 */

#ifndef __PARTICLE_H_INCLUDED__
#define __PARTICLE_H_INCLUDED__
#include <string>
#include <array>
#include <iostream>
class Particle
{
	public:
		std::string name; // an individual name
		std::string type; // determines potentials
		std::array<double, 3> position; // current position
		unsigned long int numberOfTimestep; // current time
		unsigned int count; // how many timesteps can still be skipped
		unsigned int skip; // how many timesteps in total were skipped. Determines the distribution from which new position is drawn
		double radius;	// size of the particle
		double diffusionConstant; // an individual diffusionConstant, typically determined by the type
		std::array<double, 3> cumulativeForce;	// the cumulative force for the current timestep. Should be zero after propagation

		Particle();
		~Particle();
		void move(std::array<double, 3> deviation);
		void addForce(std::array<double, 3> force);
		void resetForce();
};

#endif
