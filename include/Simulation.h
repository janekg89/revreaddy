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

class Simulation
{
	public:
		unsigned long int maxTime;
		double timestep;
		double temperature;
		double kBoltzmann;
		void run();
		void propagate(ActiveParticles * activeParticles, Random * random);
};

#endif
