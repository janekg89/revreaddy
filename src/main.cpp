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
	sim->squaredDistances.resize(sim->maxTime);

	sim->run(ap, random);

	printVector(sim->squaredDistances);

	ap->container.clear();
	delete p1;
	delete random;
	return 0;
}


/*
int main()
{
	std::array<double, 3> arr  = {1., 2., 3.};
	printArray<arr.size()>(arr);
	return 0;
}
*/
