/* main.cpp
 *
 * This is a pure C++ test environment for
 * the revreaddy python extension.
 *
 */

#include "Random.h"
#include "Particle.h"
#include "Simulation.h"
#include "SingleParticleDiffusion.h"
#include "Filehandler.h"
#include "utils.h"


int main()
{
	SingleParticleDiffusion * sim = new SingleParticleDiffusion;
	Particle * p1 = new Particle();
	const std::array<double, 3> x0 = {0., 0., 0.};
	sim->addParticle(p1);
	sim->activeParticles[0].position = x0;
	sim->initialPosition = x0;

	sim->kBoltzmann = 1.;
	sim->maxTime	= 100000;
	sim->temperature= 1.;
	sim->timestep	= 0.001;
	sim->isPeriodic = true;
	sim->boxSize = 10.;
	sim->squaredDistances.resize(sim->maxTime);
	sim->trajectory.resize(sim->maxTime);

	sim->run();

	//printVector(sim->squaredDistances);
	Filehandler * filehandler;
	filehandler->writeSingleParticleTrajectoryXyz(sim->trajectory);
	return 0;
}
