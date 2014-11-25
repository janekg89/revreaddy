/* Simulation.cpp
 * author: Christoph Froehner
 */

#include "Simulation.h"
#include "utils.h"

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
		// particleForces()
		// groupForces()
		propagate();
		this->recordObservables(t);
	}
	std::cout << "Simulation has finished\n";
}

void Simulation::propagate()
{
	std::array<double, 3> noiseTerm;
	double noisePrefactor;
	std::array<double, 3> forceTerm;
	double forcePrefactor;
	//for (int i=0; i<activeParticles.size(); i++)
	for (auto&& particle : activeParticles)
	{
		noiseTerm = random->normal3D();
		noisePrefactor = sqrt(2. * particle.diffusionConstant * timestep);
		noiseTerm[0] *= noisePrefactor;
		noiseTerm[1] *= noisePrefactor;
		noiseTerm[2] *= noisePrefactor;

		forcePrefactor = timestep * particle.diffusionConstant / (kBoltzmann * temperature);
		forceTerm[0] = particle.cumulativeForce[0] * forcePrefactor;
		forceTerm[1] = particle.cumulativeForce[1] * forcePrefactor;
		forceTerm[2] = particle.cumulativeForce[2] * forcePrefactor;

		particle.move(noiseTerm);
		particle.resetForce();

		if (isPeriodic)
		{
			if (particle.position[0] < (-0.5 * boxSize) ) {particle.position[0] += boxSize;}
			if (particle.position[0] >= (0.5 * boxSize) ) {particle.position[0] -= boxSize;}
			if (particle.position[1] < (-0.5 * boxSize) ) {particle.position[1] += boxSize;}
			if (particle.position[1] >= (0.5 * boxSize) ) {particle.position[1] -= boxSize;}
			if (particle.position[2] < (-0.5 * boxSize) ) {particle.position[2] += boxSize;}
			if (particle.position[2] >= (0.5 * boxSize) ) {particle.position[2] -= boxSize;}
		}
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
