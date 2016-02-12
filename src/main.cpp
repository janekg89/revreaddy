#include "Simulation.h"
#include "Config.h"
#define print(x) std::cout << x << std::endl;

int main()
{
	Simulation * sim = new Simulation();

	sim->config->kT= 2.437; // kJ/mol
	sim->config->timestep = 1e-10; // seconds
	sim->config->isPeriodic = false;
	sim->config->boxsize = 110.; //nanometer
	sim->isReversibleDynamics = false;
	sim->isReversibleReactions = false;
	sim->useNeighborList = true;

	sim->config->new_Type("A", 1.5 , 143e6); // diff units are nanometer^2 / second
	sim->config->new_Type("B", 3. , 71e6);

	std::vector<unsigned> affected = {0,1};
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
	print(sim->world->getNumberOfParticles())
	sim->run(100);
	return 0;
}