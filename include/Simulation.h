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
		Random * random; // the random number generator
		Potential * potential; // the force/energy handler
		std::vector<Particle> activeParticles;
		//std::vector<Particle> consideredParticles; // for later adaptive timestepping methods
		unsigned long int maxTime;
		double timestep;
		double temperature;
		double kBoltzmann;
		bool isPeriodic;
		double boxSize;

		void addParticle(std::array<double,3> initPos, double rad, double diffConst);
		void run();
		void propagate();
		void recordObservables(unsigned long int t);
		// should double loop (i,j) over activeParticles and call according Forcetype for
		// particle pair (i,j)
		void calculateRepulsionForces();

		Simulation();
		virtual ~Simulation();

};

#endif
