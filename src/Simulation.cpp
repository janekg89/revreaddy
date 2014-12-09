/* Simulation.cpp
 * author: Christoph Froehner
 */

#include "Simulation.h"
#include "utils.h"

Simulation::Simulation()
{
	random            = new Random("ranlxs0");
	potential         = new Potential();
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
	delete potential;
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
	std::array<double, 3> r_ij; //connecting vector from particle i to j
	double rSquared; //distance of particles i,j squared
	double cutoffSquared; //cutoff distance of particles i,j softcore interaction
	for (int i=0; i<activeParticles.size(); i++)
	{
		for (int j=i+1; j<activeParticles.size(); j++)
		{
			r_ij = getMinDistance(activeParticles[i].position, activeParticles[j].position);
			rSquared = r_ij[0]*r_ij[0] + r_ij[1]*r_ij[1] + r_ij[2]*r_ij[2];
			cutoffSquared = pow(activeParticles[i].radius + activeParticles[j].radius, 2.);
			forceI = potential->LJ1206(r_ij, rSquared, cutoffSquared, repulsionStrength);
			forceJ[0] = -1. * forceI[0];
			forceJ[1] = -1. * forceI[1];
			forceJ[2] = -1. * forceI[2];
			activeParticles[i].addForce(forceI);
			activeParticles[j].addForce(forceJ);
		}
	}	
}


std::array<double,3> Simulation::getMinDistance( std::array<double,3> r_i, std::array<double,3> r_j)
// Return the minimum distance vector, pointing from r_i to r_j
{
	double dx, dy, dz;
	std::array<double,3> r_ij;
	dx = r_j[0] - r_i[0];
	dy = r_j[1] - r_i[1];
	dz = r_j[2] - r_i[2];
	if ( isPeriodic )
	{
		// copysign returns a value with the magnitude of the first arg and
		// the sign of the second arg
		if ( fabs(dx) > 0.5*this->boxsize ) { dx -= copysign(this->boxsize, dx); }
 		if ( fabs(dy) > 0.5*this->boxsize ) { dy -= copysign(this->boxsize, dy); }
		if ( fabs(dz) > 0.5*this->boxsize ) { dz -= copysign(this->boxsize, dz); }
	}
	r_ij[0] = dx;
	r_ij[1] = dy;
	r_ij[2] = dz;
	return r_ij;
}

void Simulation::addParticle(std::array<double,3> initPos, double rad, double diffConst)
{
	Particle * particle = new Particle();
	particle->position = initPos;
	particle->radius = rad;
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
