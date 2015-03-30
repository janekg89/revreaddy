/* main.cpp
 *
 * This is a pure C++ test environment for
 * the revreaddy python/cython extension.
 *
 */

#include "Simulation.h"
#include "utils.h"
#include "Trajectory.h"
#include "RadialDistribution.h"

int main()
{
	Simulation * sim = new Simulation(true);

	sim->kBoltzmann = 1.;
	sim->maxTime	= 1000;
	sim->temperature= 1.;
	sim->timestep	= 0.01;
	sim->isPeriodic = true;
	sim->boxsize = 10.;
	sim->repulsionStrength = 1.;

	std::vector<double> x0 = {0., 0., 0.};
	for (double i=0; i<5; i++)
		for (double j=0; j<5; j++)
			for (double k=0; k<5; k++) {
				x0[0] = -0.49 * sim->boxsize + i;
				x0[1] = -0.49 * sim->boxsize + j;
				x0[2] = -0.49 * sim->boxsize + k;
				sim->addParticle(x0, 2);
			}
	sim->new_Trajectory("traj.xyz");
	sim->new_MeanSquaredDisplacement("msd.dat", 2);
	sim->run();
	sim->writeAllObservablesToFile();
	sim->deleteAllObservables();
	return 0;
}
