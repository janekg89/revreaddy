/* Simulation.cpp
 * author: Christoph Froehner
 */

#include "Simulation.h"
#include "utils.h"

Simulation::Simulation()
{
	random = new Random("ranlxs0");
	potential = new Potential();
}

Simulation::~Simulation()
{
	delete random;
	delete potential;
}

void Simulation::run()
{
	std::cout << "Simulation has started\n";
	this->recordObservables(0);
	for (unsigned long int t = 1; t < maxTime; t++)
	{
		// particleForces()
		// groupForces()
		calculateForces();
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

void Simulation::calculateForces()
{
	std::array<double, 3> forceI;
	std::array<double, 3> forceJ;
	std::array<double, 3> r_ij; //connecting vector from particle i to j
	double rSquared; //distance of particles i,j squared
	double cutoffSquared; //cutoff distance of particles i,j softcore interaction
	double strength = 2.; // force constant for the softcore interaction
	for (int i=0; i<activeParticles.size(); i++)
	{
		for (int j=i; j<activeParticles.size(); j++)
		{
			r_ij[0] = activeParticles[j].position[0] - activeParticles[i].position[0];
			r_ij[1] = activeParticles[j].position[1] - activeParticles[i].position[1];
			r_ij[2] = activeParticles[j].position[2] - activeParticles[i].position[2];
			rSquared = r_ij[0]*r_ij[0] + r_ij[1]*r_ij[1] + r_ij[2]*r_ij[2];
			cutoffSquared = pow(activeParticles[i].radius + activeParticles[j].radius, 2.);
			if (rSquared < cutoffSquared)
			{
				forceI = potential->softcoreForce(r_ij, rSquared, cutoffSquared, strength);
				forceJ[0] = -1. * forceI[0];
				forceJ[1] = -1. * forceI[1];
				forceJ[2] = -1. * forceI[2];
				activeParticles[i].addForce(forceI);
				activeParticles[j].addForce(forceJ);
			}
		}
	}	
}
void Simulation::addParticle(Particle * particle)
{
	activeParticles.push_back(*particle);
}

void Simulation::recordObservables(unsigned long int)
{
	//std::cout << "recordObservables of Simulation called\n";
}
