/* main.cpp
 *
 * This is a pure C++ test environment for
 * the revreaddy python extension.
 *
 */

#include "Random.h"
#include "Particle.h"
#include "ActiveParticles.h"
#include "Simulation.h"
#include "SingleParticleDiffusion.h"
#include "utils.h"

int main()
{
	Particle * p1 = new Particle();
	Random * random = new Random("ranlxs0");
	ActiveParticles * ap = new ActiveParticles();
	ap->addParticle(p1);

	SingleParticleDiffusion * sim = new SingleParticleDiffusion;
	sim->kBoltzmann = 1.;
	sim->maxTime	= 1000;
	sim->temperature= 1.;
	sim->timestep	= 0.001;
	sim->meanSquaredDistances = new double[sim->maxTime];

	sim->run(ap, random);

	printArray(sim->meanSquaredDistances, sim->maxTime);

	ap->container.clear();
	delete p1;
	delete random;
	return 0;
}
