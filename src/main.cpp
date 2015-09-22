/* main.cpp
 *
 * This is a pure C++ test environment for
 * the revreaddy python/cython extension.
 *
 */

#include "Simulation.h"
#include "Config.h"

int main()
{
	Simulation * sim = new Simulation();

	sim->config->kBoltzmann = 1.;
	sim->config->maxTime	= 1000;
	sim->config->temperature= 1.;
	sim->config->timestep	= 0.001;
	sim->config->isPeriodic = true;
	sim->config->boxsize = 10.;
	sim->config->isReversibleDynamics = true;
	sim->config->isReversibleReactions = true;

	sim->config->new_Type("softy", 0.5 , 1., 0.5);
	std::vector<unsigned int> affected = {0, 0};
	sim->config->new_SoftRepulsion("softy repulsion", affected, 1.);
	
	std::vector<double> x0 = {0., 0., 0.};
	for (double i=0; i<5; i++)
		for (double j=0; j<5; j++)
			for (double k=0; k<5; k++) {
				x0[0] = -0.49 * sim->config->boxsize + i;
				x0[1] = -0.49 * sim->config->boxsize + j;
				x0[2] = -0.49 * sim->config->boxsize + k;
				sim->world->addParticle(x0, 0);
			}

	sim->config->new_Trajectory(1, "traj.xyz");
	sim->config->new_MeanSquaredDisplacement(1, "msd.dat", 0);
	sim->run();
	sim->config->writeAllObservablesToFile();
	sim->config->deleteAllObservables();

	std::cout << "acc dynamics " << sim->world->acceptionsDynamics << std::endl;
	std::cout << "rej dynamics " << sim->world->rejectionsDynamics << std::endl;
	return 0;
}