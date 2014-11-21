/* Simulation.h
 * author: Christoph Froehner
 *
 * Simulation performs the time loop which propagates
 * the particles.
 */

#ifndef __SIMULATION_H_INCLUDED__
#define __SIMULATION_H_INCLUDED__
#include <math.h>
#include "Particle.h"
#include "ActiveParticles.h"
#include "Random.h"
#include <array>
#include <vector>

// TODO extend this class to omit ActiveParticles class

class Simulation
{
	public:
		unsigned long int maxTime;
		double timestep;
		double temperature;
		std::vector<Particle> activeParticles;
		//std::vector<Particle> consideredParticles; // for later adaptive timestepping methods
		void run();
		void timeloop();
		void propagate(ActiveParticles * activeParticles, Random * random);
		// should double loop (i,j) over activeParticles and call according Forcetype for
		// particle pair (i,j)
		void calculateForces();

};

#endif
