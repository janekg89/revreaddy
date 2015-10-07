#include "Simulation.h"
#include "Config.h"
#define print(x) std::cout << x << std::endl;

int main()
{
	Simulation * sim = new Simulation();

	sim->config->kBoltzmann = 1.;
	sim->config->temperature= 2.437;
	sim->config->timestep	= 1e-10;
	sim->config->isPeriodic = false;
	sim->config->boxsize = 110.;
	sim->config->isReversibleDynamics = false;
	sim->config->isReversibleReactions = false;
	sim->config->maxTime	= 5;

	sim->config->useNeighborList = true;
	sim->config->numberBoxes = 36;

	sim->config->new_Type("A", 1.5 , 143e6, 1.5);
	sim->config->new_Type("B", 3. , 71e6, 3.);

	std::vector<unsigned> affected = {0,1,2};
	double kappa = 5.;
	sim->config->deleteAllGeometries();
	sim->config->new_Wall(std::vector<double> {1.,0.,0.}, std::vector<double> {-50.,0.,0.}, kappa, affected);
	sim->config->new_Wall(std::vector<double> {-1.,0.,0.}, std::vector<double> {50.,0.,0.}, kappa, affected);
	sim->config->new_Wall(std::vector<double> {0.,1.,0.}, std::vector<double> {0.,-50.,0.}, kappa, affected);
	sim->config->new_Wall(std::vector<double> {0.,-1.,0.}, std::vector<double> {0.,50.,0.}, kappa, affected);
	sim->config->new_Wall(std::vector<double> {0.,0.,1.}, std::vector<double> {0.,0.,-50.}, kappa, affected);
	sim->config->new_Wall(std::vector<double> {0.,0.,-1.}, std::vector<double> {0.,0.,50.}, kappa, affected);	
	
	std::vector<double> x0 = {0., 0., 0.};
	for (double i=0; i<95; i+=7.063305534550977)
		for (double j=0; j<95; j+=7.063305534550977)
			for (double k=0; k<95; k+=7.063305534550977) {
				x0[0] = -47. + i;
				x0[1] = -47. + j;
				x0[2] = -47. + k;
				sim->world->addParticle(x0, 0);
			}
	print(sim->config->getParticleNumber())
	sim->run();
	return 0;
}