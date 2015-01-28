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
	energy            = 0.;
	oldEnergy         = 0.;
	acceptions        = 0;
	rejections        = 0;
}

Simulation::~Simulation()
{
	delete random;
	delete force;
}

void Simulation::run()
{
	std::cout << "Simulation has started\n";
	this->resetForces();
	this->calculateRepulsionForcesEnergies();
	this->recordObservables(0);
	for (unsigned long int t = 1; t < maxTime; t++)
	{
		// groupForces()
		this->saveOldState();
		this->propagate(); // propose
		this->resetForces();
		this->calculateRepulsionForcesEnergies(); // calculate energy and force
		this->acceptOrReject();
		this->recordObservables(t);
	}
	std::cout << "Simulation has finished\n";
}

void Simulation::saveOldState()
{
	this->oldEnergy = this->energy;
	for (int i=0; i<activeParticles.size(); i++) {
		activeParticles[i].oldPosition = activeParticles[i].position;
		activeParticles[i].oldForce    = activeParticles[i].cumulativeForce;
	}
}

void Simulation::propagate()
{
	std::vector<double> noiseTerm = {0.,0.,0.};
	std::vector<double> forceTerm = {0.,0.,0.};
	double noisePrefactor;
	double forcePrefactor;
	for (int i=0; i<activeParticles.size(); i++)
	{
		noiseTerm = random->normal3D();
		noisePrefactor = sqrt(
			2. * activeParticles[i].diffusionConstant * timestep);
		noiseTerm[0] *= noisePrefactor;
		noiseTerm[1] *= noisePrefactor;
		noiseTerm[2] *= noisePrefactor;

		forcePrefactor = timestep * activeParticles[i].diffusionConstant 
			/ (kBoltzmann * temperature);
		forceTerm[0] = activeParticles[i].cumulativeForce[0] * forcePrefactor;
		forceTerm[1] = activeParticles[i].cumulativeForce[1] * forcePrefactor;
		forceTerm[2] = activeParticles[i].cumulativeForce[2] * forcePrefactor;

		activeParticles[i].move(noiseTerm);
		activeParticles[i].move(forceTerm);

		if (isPeriodic)
		{
			if (activeParticles[i].position[0] < (-0.5 * boxsize) ) {
				activeParticles[i].position[0] += boxsize;}
			else if (activeParticles[i].position[0] >= (0.5 * boxsize) ) {
				activeParticles[i].position[0] -= boxsize;}
			if (activeParticles[i].position[1] < (-0.5 * boxsize) ) {
				activeParticles[i].position[1] += boxsize;}
			else if (activeParticles[i].position[1] >= (0.5 * boxsize) ) {
				activeParticles[i].position[1] -= boxsize;}
			if (activeParticles[i].position[2] < (-0.5 * boxsize) ) {
				activeParticles[i].position[2] += boxsize;}
			else if (activeParticles[i].position[2] >= (0.5 * boxsize) ) {
				activeParticles[i].position[2] -= boxsize;}
		}
	}
}

void Simulation::recordObservables(unsigned long int t)
{
	for (auto* obs : this->observables) {
		obs->record(this->activeParticles, t);
	}
}

/* The reversible version of this function also calculates energies.*/
void Simulation::calculateRepulsionForcesEnergies()
{
	std::vector<double> forceI = {0.,0.,0.};
	std::vector<double> forceJ = {0.,0.,0.};
	// connecting vector from particle i to j
	std::vector<double> r_ij   = {0.,0.,0.}; 
	double rSquared     = 0.8; // distance of particles i,j squared
	double radiiSquared = 1.; // squared sum of particles i,j radii
	double energyBuffer = 0.; // interaction energy of particle pair (i,j)
	this->energy = 0.;
	for (int i=0; i<activeParticles.size(); i++) {
		for (int j=i+1; j<activeParticles.size(); j++) {
			getMinDistanceVector(
				r_ij,
				activeParticles[i].position, 
				activeParticles[j].position, 
				this->isPeriodic, 
				this->boxsize);
			rSquared = r_ij[0]*r_ij[0] + r_ij[1]*r_ij[1] + r_ij[2]*r_ij[2];
			radiiSquared = pow(
				activeParticles[i].radius + activeParticles[j].radius,
				2.);
			force->repulsionForceEnergy(
				forceI,
				energyBuffer,
				r_ij,
				rSquared,
				radiiSquared,
				repulsionStrength, 
				activeParticles[i].type,
				activeParticles[j].type);
			forceJ[0] = -1. * forceI[0];
			forceJ[1] = -1. * forceI[1];
			forceJ[2] = -1. * forceI[2];
			activeParticles[i].addForce(forceI);
			activeParticles[j].addForce(forceJ);
			this->energy += energyBuffer;
		}
	}
}

void Simulation::resetForces()
{
	for (int i=0; i<activeParticles.size(); i++) {
		activeParticles[i].resetForce();
	}
}

void Simulation::acceptOrReject()
{
	double acceptance = 1.;
	double firstTerm  = 0.;
	double secondTerm = 0.;
	std::vector<double> deltaX = {0.,0.,0.};
	for (int i=0; i<activeParticles.size(); i++) {
		getMinDistanceVector(
			deltaX,
			activeParticles[i].oldPosition,
			activeParticles[i].position,
			this->isPeriodic,
			this->boxsize);
		firstTerm  += deltaX[0]
		            * ( activeParticles[i].oldForce[0]
		              + activeParticles[i].cumulativeForce[0] );
		firstTerm  += deltaX[1]
		            * ( activeParticles[i].oldForce[1]
		              + activeParticles[i].cumulativeForce[1] );
		firstTerm  += deltaX[2]
		            * ( activeParticles[i].oldForce[2]
		              + activeParticles[i].cumulativeForce[2] );
		secondTerm += activeParticles[i].diffusionConstant
		            * ( activeParticles[i].cumulativeForce[0]
		              * activeParticles[i].cumulativeForce[0]
		              + activeParticles[i].cumulativeForce[1]
		              * activeParticles[i].cumulativeForce[1]
		              + activeParticles[i].cumulativeForce[2]
		              * activeParticles[i].cumulativeForce[2]
		              + activeParticles[i].oldForce[0]
		              * activeParticles[i].oldForce[0]
		              + activeParticles[i].oldForce[1]
		              * activeParticles[i].oldForce[1]
		              + activeParticles[i].oldForce[2]
		              * activeParticles[i].oldForce[2] );
	}
	firstTerm  *= 0.5;
	secondTerm *= this->timestep / (4. * this->kBoltzmann * this->temperature);
	acceptance = firstTerm + secondTerm + this->energy - this->oldEnergy;
	acceptance /= -1. * this->kBoltzmann * this->temperature;
	acceptance = exp( acceptance );
	// accept or reject
	// std::cout << "acceptance "<<acceptance<<"\n";
	if ( acceptance > 1. ) {
		/*accept = do nothing. particles keep their new positions and forces*/
		acceptions += 1;
	}
	else {
		double uniform = this->random->uniform();
		if ( uniform < acceptance ) {/* accept */ acceptions += 1;}
		else {/* reject = restore old positions, forces, energy */
			rejections += 1;
			this->energy = this->oldEnergy;
			for (int i=0; i<activeParticles.size(); i++) {
				activeParticles[i].position 
					= activeParticles[i].oldPosition;
				activeParticles[i].cumulativeForce
					= activeParticles[i].oldForce;
			}
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
