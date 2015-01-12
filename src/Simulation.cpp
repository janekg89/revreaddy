/* Simulation.cpp
 * author: Christoph Froehner
 */

#include "Simulation.h"

Simulation::Simulation()
{
	random            = new Random("ranlxs0");
	force             = new Force();
	timestep          = 0.001;
	temperature       = 1.;
	kBoltzmann        = 1.;
	repulsionStrength = 1.;
	isPeriodic        = true;
	boxsize           = 10.;
}

Simulation::~Simulation()
{
	delete random;
	delete force;
}

void Simulation::run()
{
	std::cout << "Simulation has started\n";
	this->recordObservables(0);
	for (unsigned long int t = 1; t < maxTime; t++)
	{
		// groupForces()
		this->calculateRepulsionForces();
		this->propagate();
		this->recordObservables(t);
	}
	std::cout << "Simulation has finished\n";
}

void Simulation::propagate()
{
	std::array<double, 3> noiseTerm;
	std::array<double, 3> forceTerm;
	double noisePrefactor;
	double forcePrefactor;
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
		particle.move(forceTerm);
		particle.resetForce();

		if (isPeriodic)
		{
			if (particle.position[0] < (-0.5 * boxsize) ) {particle.position[0] += boxsize;}
			if (particle.position[0] >= (0.5 * boxsize) ) {particle.position[0] -= boxsize;}
			if (particle.position[1] < (-0.5 * boxsize) ) {particle.position[1] += boxsize;}
			if (particle.position[1] >= (0.5 * boxsize) ) {particle.position[1] -= boxsize;}
			if (particle.position[2] < (-0.5 * boxsize) ) {particle.position[2] += boxsize;}
			if (particle.position[2] >= (0.5 * boxsize) ) {particle.position[2] -= boxsize;}
		}
	}
}

void Simulation::calculateRepulsionForces()
{
	std::array<double, 3> forceI;
	std::array<double, 3> forceJ;
	std::array<double, 3> r_ij; // connecting vector from particle i to j
	double rSquared; // distance of particles i,j squared
	double radiiSquared; // squared sum of particles i,j radii
	for (int i=0; i<activeParticles.size(); i++)
	{
		for (int j=i+1; j<activeParticles.size(); j++)
		{
			r_ij = getMinDistanceVector(
				activeParticles[i].position, 
				activeParticles[j].position, 
				this->isPeriodic, 
				this->boxsize
			);
			rSquared = r_ij[0]*r_ij[0] + r_ij[1]*r_ij[1] + r_ij[2]*r_ij[2];
			radiiSquared = pow(activeParticles[i].radius + activeParticles[j].radius, 2.);
			try {
				forceI = force->repulsion(r_ij, rSquared, radiiSquared, repulsionStrength, 
					activeParticles[i].type, activeParticles[j].type);
			}
			catch(const char* msg) {
				std::cerr << msg << "\n";
				forceI = {0.,0.,0.};
			}
			forceJ[0] = -1. * forceI[0];
			forceJ[1] = -1. * forceI[1];
			forceJ[2] = -1. * forceI[2];
			activeParticles[i].addForce(forceI);
			activeParticles[j].addForce(forceJ);
		}
	}	
}

void Simulation::addParticle(std::vector<double> initPos, std::string particleType, double rad, double diffConst)
{
	// first convert initPos from vector to array
	std::array<double,3> initPosArr;
	try {
		if ( initPos.size() == 3 ) {
			initPosArr[0] = initPos[0];
			initPosArr[1] = initPos[1];
			initPosArr[2] = initPos[2];
		}
		else {
			throw   "Particles' initial position has dimension mismatch! \
					Particle will be placed at {0,0,0}.";	
		}
	}
	catch(const char* msg) {
		std::cerr << msg << "\n";
		initPosArr = {0., 0., 0.};
	}
	Particle * particle         = new Particle();
	particle->position          = initPosArr;
	particle->type              = particleType;
	particle->radius            = rad;
	particle->diffusionConstant = diffConst;
	this->activeParticles.push_back(*particle);//push_back copies arg into vec
	delete particle;
}

void Simulation::recordObservables(unsigned long int t)
{
	for (auto* obs : observables)
	{
		obs->record(this->activeParticles, t);
	}
}
