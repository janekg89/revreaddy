/* Simulation.cpp
 * author: Christoph Froehner
 */

#include "Simulation.h"

void Simulation::propagate(ActiveParticles * activeParticles, Random * random)
{
	double noiseTerm[3];
	double noisePrefactor;
	double forceTerm[3];
	double forcePrefactor;
	for (int i=0; i<activeParticles->container.size(); i++)
	{
		random->normal3D(noiseTerm);
		noisePrefactor = sqrt(2. * activeParticles->container[i].diffusionConstant * this->timestep);
		noiseTerm[0] *= noisePrefactor;
		noiseTerm[1] *= noisePrefactor;
		noiseTerm[2] *= noisePrefactor;

		forcePrefactor = this->timestep * activeParticles->container[i].diffusionConstant / (this->kBoltzmann * this->temperature);
		forceTerm[0] = activeParticles->container[i].force[0] * forcePrefactor;
		forceTerm[1] = activeParticles->container[i].force[1] * forcePrefactor;
		forceTerm[2] = activeParticles->container[i].force[2] * forcePrefactor;

		activeParticles->container[i].move(noiseTerm);
		activeParticles->container[i].resetForce();
	}
}
