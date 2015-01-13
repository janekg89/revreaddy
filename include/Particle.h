/* Particle.h
 * author: Christoph Froehner
 *
 * This data structure will hold information
 * about a single particle. It will be constructed
 * for every particle in the simulation.
 */

#ifndef __PARTICLE_H_INCLUDED__
#define __PARTICLE_H_INCLUDED__
#include <string>
#include <iostream>
#include <vector>
class Particle
{
	public:
		std::string type; // determines potentials
		std::vector<double> position; // current position
		unsigned int count; // how many timesteps can still be skipped
		/* how many timesteps in total were skipped. Determines 
		 * the distribution from which new position is drawn */
		unsigned int skip;
		double radius;	// size of the particle
		/* an individual diffusionConstant, typically determined by the type */
		double diffusionConstant; 
		/* the cumulative force for the current timestep. 
		 * Should be zero after propagation */
		std::vector<double> cumulativeForce;

		Particle();
		~Particle();
		void move(std::vector<double> deviation);
		void addForce(std::vector<double> force);
		void resetForce();
};

#endif
