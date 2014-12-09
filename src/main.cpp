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
	Simulation * sim = new Simulation();

	sim->kBoltzmann = 1.;
	sim->maxTime	= 10000;
	sim->temperature= 1.;
	sim->timestep	= 0.001;
	sim->isPeriodic = true;
	sim->boxsize = 8.;
	sim->repulsionStrength = 2000.;

	std::array<double,3> x0 = {0.,0.,0.};
	for (int i=0; i<100; i++) {
		sim->addParticle(x0, 1., 1.0);
	}
	
	sim->run();

	Trajectory * traj = new Trajectory();
	RadialDistribution * rad = new RadialDistribution( 100 , sim );
	std::vector<double> ranges;
	for (double i=0; i<101; i++) {
		ranges.push_back(0.1*i);
	}
	rad->setRange(ranges);
	sim->observables.push_back(traj);
	sim->observables.push_back(rad);

	sim->run();

	traj->writeBufferToFile();
	rad->writeBufferToFile();
	return 0;
}
