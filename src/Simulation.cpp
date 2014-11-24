/* Simulation.cpp
 * author: Christoph Froehner
 */

#include "Simulation.h"

Simulation::Simulation()
{
	random = new Random("ranlxs0");
}

Simulation::~Simulation()
{
	delete random;
}


void Simulation::run()
{
	std::cout << "Simulation has started\n";
	this->recordObservables(0);
	for (unsigned long int t = 1; t < maxTime; t++)
	{
		propagate();
		this->recordObservables(t);
	}
}

void Simulation::propagate()
{
	std::array<double, 3> noiseTerm;
	double noisePrefactor;
	std::array<double, 3> forceTerm;
	double forcePrefactor;
	for (int i=0; i<activeParticles.size(); i++)
	{
		noiseTerm = random->normal3D();
		noisePrefactor = sqrt(2. * activeParticles[i].diffusionConstant * timestep);
		noiseTerm[0] *= noisePrefactor;
		noiseTerm[1] *= noisePrefactor;
		noiseTerm[2] *= noisePrefactor;

		forcePrefactor = timestep * activeParticles[i].diffusionConstant / (kBoltzmann * temperature);
		forceTerm[0] = activeParticles[i].cumulativeForce[0] * forcePrefactor;
		forceTerm[1] = activeParticles[i].cumulativeForce[1] * forcePrefactor;
		forceTerm[2] = activeParticles[i].cumulativeForce[2] * forcePrefactor;

		activeParticles[i].move(noiseTerm);
		activeParticles[i].resetForce();
	}
}

void Simulation::addParticle(Particle * particle)
{
	activeParticles.push_back(*particle);
}


void Simulation::recordObservables(unsigned long int)
{
	std::cout << "recordObservables of Simulation called\n";
}

