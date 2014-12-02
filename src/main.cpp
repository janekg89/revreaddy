/* main.cpp
 *
 * This is a pure C++ test environment for
 * the revreaddy python/cython extension.
 *
 */

#include "Simulation.h"
#include "Filehandler.h"
#include "utils.h"
#include "Trajectory.h"

int main()
{
	Simulation * sim = new Simulation();

	sim->kBoltzmann = 1.;
	sim->maxTime	= 1000;
	sim->temperature= 1.;
	sim->timestep	= 0.001;
	sim->isPeriodic = true;
	sim->boxsize = 20.;
	sim->repulsionStrength = 200.;

	std::array<double,3> x0 = {0.,0.,0.};
	for (int i=0; i<30; i++) {
		sim->addParticle(x0, 4., 1.);
	}
	
	Trajectory * traj = new Trajectory();

	sim->observables.push_back(traj);
	sim->run();

	traj->writeBufferToFile();

	printArray(sim->activeParticles[0].position);
	printArray(sim->activeParticles[1].position);
	return 0;
}
