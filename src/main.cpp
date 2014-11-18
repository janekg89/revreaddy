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
	Particle * p1 = new Particle();
	Random * random = new Random("ranlxs0");
	print3DArray(p1->position);

	std::array<double, 3> r = random->normal3D();
	print3DArray(r);
	p1->move(r);
	print3DArray(p1->position);

	ActiveParticles * ap = new ActiveParticles();
	ap->addParticle(p1);

	print3DArray(ap->container[0].position);

	return 0;
}
*/
