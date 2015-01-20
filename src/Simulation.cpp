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
	verbose           = false;
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
	std::vector<double> noiseTerm = {0.,0.,0.};
	std::vector<double> forceTerm = {0.,0.,0.};
	double noisePrefactor;
	double forcePrefactor;
	for (auto&& particle : activeParticles)
	{
		noiseTerm = random->normal3D();
		noisePrefactor = sqrt(2. * particle.diffusionConstant * timestep);
		noiseTerm[0] *= noisePrefactor;
		noiseTerm[1] *= noisePrefactor;
		noiseTerm[2] *= noisePrefactor;

		forcePrefactor = timestep * particle.diffusionConstant 
			/ (kBoltzmann * temperature);
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

// TODO avoid alocation in every function call. 
// Make new Simulation variables instead.
void Simulation::calculateRepulsionForces()
{
	std::vector<double> forceI = {0.,0.,0.};
	std::vector<double> forceJ = {0.,0.,0.};
	// connecting vector from particle i to j
	std::vector<double> r_ij = {0.,0.,0.}; 
	double rSquared; // distance of particles i,j squared
	double radiiSquared; // squared sum of particles i,j radii
	for (int i=0; i<activeParticles.size(); i++) {
		for (int j=i+1; j<activeParticles.size(); j++) {
			r_ij = getMinDistanceVector(
				activeParticles[i].position, 
				activeParticles[j].position, 
				this->isPeriodic, 
				this->boxsize);
			rSquared = r_ij[0]*r_ij[0] + r_ij[1]*r_ij[1] + r_ij[2]*r_ij[2];
			radiiSquared = pow(
				activeParticles[i].radius + activeParticles[j].radius,
				2.);
			try {
				forceI = force->repulsionForce(
					r_ij,
					rSquared,
					radiiSquared,
					repulsionStrength, 
					activeParticles[i].type,
					activeParticles[j].type);
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

void Simulation::addParticle(
	std::vector<double> initPos,
	std::string particleType,
	double rad,
	double diffConst)
{
	Particle * particle         = new Particle();
	try {
		if ( initPos.size() == 3 ) { particle->position  = initPos; }
		else {
			throw "Particles' initial position has dimension mismatch!" 
			      "Particle will be placed at {0,0,0}.";	
		}
	}
	catch(const char* msg) {
		std::cerr << msg << "\n";
		particle->position = {0., 0., 0.};
	}
	particle->type              = particleType;
	particle->radius            = rad;
	particle->diffusionConstant = diffConst;
	this->activeParticles.push_back(*particle);//push_back copies arg into vec
	delete particle;
	if (this->verbose) {std::cout << "Particle added.\n";}
}

void Simulation::recordObservables(unsigned long int t)
{
	for (auto* obs : this->observables)
	{
		obs->record(this->activeParticles, t);
	}
}

std::vector<double> Simulation::getPosition(int index)
{
	return this->activeParticles[index].position;
}

void Simulation::setPosition(int index, std::vector<double> newPos)
{
	if (newPos.size() == 3) {
		this->activeParticles[index].position[0] = newPos[0];
		this->activeParticles[index].position[1] = newPos[1];
		this->activeParticles[index].position[2] = newPos[2];
	}
	else {
		throw "New position has dimension mismatch!"
		      "Particle remains at its old position.";
	}
}

int Simulation::getParticleNumber()
{
	return this->activeParticles.size();
}

void Simulation::deleteAllParticles()
{
	 this->activeParticles.erase(
		this->activeParticles.begin(),
		this->activeParticles.begin() + this->activeParticles.size()
	);
}

void Simulation::writeAllObservablesToFile()
{
	for (auto* obs : this->observables) {
		obs->writeBufferToFile();
	}
}

std::string Simulation::showObservables()
{
	std::string content = "Observables: ";
	if (this->observables.size() > 0) {
		for (auto* obs : this->observables) {
			content += std::string (typeid(*obs).name()) + " ";
		}
	}
	else {content += "empty";}
	content += "\n";
	return content;
}

void Simulation::deleteAllObservables()
{
	/* first delete the obs, since they we're allocated with 'new'
	 * then erase the pointers from the vector */
	for (auto* obs : this->observables) {
		delete obs;
	}
	this->observables.erase(
		this->observables.begin(),
		this->observables.begin() + this->observables.size()
	);
}

void Simulation::new_Trajectory(std::string filename)
{
	Trajectory * obs = new Trajectory();
	obs->filename = filename;
	observables.push_back(obs);
}

void Simulation::new_TrajectorySingle()
{
	TrajectorySingle * obs = new TrajectorySingle();
	observables.push_back(obs);
}

std::vector< std::vector<double> > Simulation::getTrajectorySingle()
{
	/* find the first Observable of type TrajectorySingle and return it */
	for (auto* obs : this->observables) {
		if (std::string (typeid(*obs).name()) == "16TrajectorySingle") {
			return obs->trajectory;
		}
	}
	// if that fails return a default 'zero' vector.
	std::vector< std::vector<double> > zero;
	std::vector<double> x = {0.,0.,0.};
	zero.push_back(x);
	return zero;
}
