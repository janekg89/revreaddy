/* main.cpp
 *
 * This is a pure C++ test environment for
 * the revreaddy python/cython extension.
 *
 */

#include "Simulation.h"

int main()
{
	Simulation * sim = new Simulation();

	sim->kBoltzmann = 1.;
	sim->maxTime	= 1000;
	sim->temperature= 1.;
	sim->timestep	= 0.001;
	sim->isPeriodic = true;
	sim->boxsize = 10.;
	sim->isReversibleDynamics = true;
	sim->isReversibleReactions = true;

	sim->new_Type("softy", 0.5 , 1., 0.5);
	std::vector<unsigned int> affected = {0, 0};
	sim->new_SoftRepulsion("softy repulsion", affected, 1.);
	
	std::vector<double> x0 = {0., 0., 0.};
	for (double i=0; i<5; i++)
		for (double j=0; j<5; j++)
			for (double k=0; k<5; k++) {
				x0[0] = -0.49 * sim->boxsize + i;
				x0[1] = -0.49 * sim->boxsize + j;
				x0[2] = -0.49 * sim->boxsize + k;
				sim->addParticle(x0, 0);
			}
	sim->new_Trajectory(1, "traj.xyz");
	sim->new_MeanSquaredDisplacement(1, "msd.dat", 0);
	sim->run();
	sim->writeAllObservablesToFile();
	sim->deleteAllObservables();

	std::cout << "acc dynamics " << sim->acceptionsDynamics << std::endl;
	std::cout << "rej dynamics " << sim->rejectionsDynamics << std::endl;

	return 0;
}
