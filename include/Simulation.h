/* Simulation.h
 * author: Christoph Froehner
 *
 * This is the main class, which uses all other classes.
 * Simulation performs the time loop which propagates
 * the particles.
 */

#ifndef __SIMULATION_H_INCLUDED__
#define __SIMULATION_H_INCLUDED__
#include <math.h>
#include <cmath>
#include "Particle.h"
#include "Random.h"
#include "Potential.h"
#include <array>
#include <vector>
#include <iostream>

class Simulation
{
	public:
		unsigned long int maxTime;
		double timestep;
		double temperature;
		double kBoltzmann;
		std::vector<Particle> activeParticles;
		Random * random; // the random number generator
		//std::vector<Particle> consideredParticles; // for later adaptive timestepping methods
		bool isPeriodic;
		double boxSize;
		Potential * potential; // the force/energy handler

		void addParticle(Particle * particle);
		void run();
		void propagate();
		virtual void recordObservables(unsigned long int t);
		// should double loop (i,j) over activeParticles and call according Forcetype for
		// particle pair (i,j)
		void calculateForces();

		Simulation();
		virtual ~Simulation();

};

#endif
