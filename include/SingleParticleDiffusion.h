/* SingleParticleDiffusion.h
 * author: Christoph Froehner
 *
 * The aim of this simulation is to prove obediance
 * to the diffusion equation.
 */

#ifndef __SINGLEPARTICLEDIFFUSION_H_INCLUDED__
#define __SINGLEPARTICLEDIFFUSION_H_INCLUDED__
#include "Simulation.h"

class SingleParticleDiffusion: public Simulation
{
	public:
		double* meanSquaredDistances;
		void run(ActiveParticles * activeParticles, Random * random);
};

#endif
