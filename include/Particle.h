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
#include <iostream>
#include <stdio.h>
#include <string>

class Particle
{
	public:
		std::string name; // an individual name
		std::string type; // determines potentials
		double* position; // current position
		unsigned long int timestep; // current time
		unsigned int count; // how many timesteps can still be skipped
		unsigned int skip; // how many timesteps in total were skipped. Determines the distribution from which new position is drawn
		double radius;	// size of the particle
		double diffusionConstant; // an individual diffusionConstant, typically determined by the type
		double* force;	// the cumulative force for the current timestep. Should be zero after propagation

};

#endif
