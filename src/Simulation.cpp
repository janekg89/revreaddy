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
		this->calculateRepulsionForcesRev();
		this->recordObservables(t);
		this->propagateRev();
		this->updateOldPositions();
		this->resetSingleEnergies();
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

void Simulation::propagateRev()
{
	std::vector<double> noiseTerm   = {0.,0.,0.};
	std::vector<double> forceTerm   = {0.,0.,0.};
	std::vector<double> oldForce    = {0.,0.,0.};
	std::vector<double> newForce    = {0.,0.,0.};
	std::vector<double> forceBuffer = {0.,0.,0.};
	std::vector<double> r_ij        = {0.,0.,0.};
	std::vector<double> deltaX      = {0.,0.,0.};
	double noisePrefactor;
	double forcePrefactor;
	double newSingleEnergy;
	double energyBuffer;
	double distanceBuffer;
	double radiiSquared;
	double deltaEnergy;
	double acceptance;
	double uniform;
	for (int i=0; i<activeParticles.size(); i++) {
		oldForce        = activeParticles[i].cumulativeForce;

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
		activeParticles[i].resetForce();

		/* Note that actually one would have to call getMinDistanceVector
		 * because of periodic boundary conditions. But we calc deltaX 
		 * before applying those. */
		deltaX[0] = activeParticles[i].position[0]
		          - activeParticles[i].oldPosition[0];
		deltaX[1] = activeParticles[i].position[1]
		          - activeParticles[i].oldPosition[1];
		deltaX[2] = activeParticles[i].position[2]
		          - activeParticles[i].oldPosition[2];

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
		// new position has been proposed. now calculate new energy
		// and new force
		newSingleEnergy = 0.;
		newForce = {0.,0.,0.};
		for (int j=0; j<activeParticles.size(); j++) {
			if (i == j) { }
			else {
				getMinDistanceVector(
					r_ij,
					activeParticles[i].position,
					activeParticles[j].oldPosition,//use oldPosition!
					this->isPeriodic,
					this->boxsize);
				distanceBuffer=r_ij[0]*r_ij[0]+r_ij[1]*r_ij[1]+r_ij[2]*r_ij[2];
				radiiSquared = pow(
					activeParticles[i].radius + activeParticles[j].radius,
					2.);
				force->repulsionForceEnergy(
					forceBuffer,
					energyBuffer,
					r_ij,
					distanceBuffer,
					radiiSquared,
					repulsionStrength,
					activeParticles[i].type,
					activeParticles[j].type);
				newSingleEnergy += energyBuffer;
				newForce[0] += forceBuffer[0];
				newForce[1] += forceBuffer[1];
				newForce[2] += forceBuffer[2];
			}
		}
		// Calculate acceptance probability
		deltaEnergy = 0.5*(newSingleEnergy - activeParticles[i].singleEnergy);
		acceptance  = 0.5*( deltaX[0]*(newForce[0]+oldForce[0])
		                  + deltaX[1]*(newForce[1]+oldForce[1])
		                  + deltaX[2]*(newForce[2]+oldForce[2]));
		acceptance += activeParticles[i].diffusionConstant * this->timestep
			* ( newForce[0]*newForce[0] + newForce[1]*newForce[1]
			  + newForce[2]*newForce[2] + oldForce[0]*oldForce[0]
			  + oldForce[1]*oldForce[1] + oldForce[2]*oldForce[2] )
			/ (4. * this->kBoltzmann * this->temperature);
		acceptance += deltaEnergy;
		acceptance /= -1. * this->kBoltzmann * this->temperature;
		acceptance = exp( acceptance );
		// accept or reject
		if ( acceptance > 1. ) {
			/* accept = do nothing. particle is already moved */
		}
		else {
			uniform = random->uniform();
			if ( uniform < acceptance ) {/* accept */}
			else {/* reject */
				activeParticles[i].position = activeParticles[i].oldPosition;
			}
		}
	}
}

void Simulation::calculateRepulsionForces()
{
	std::vector<double> forceI = {0.,0.,0.};
	std::vector<double> forceJ = {0.,0.,0.};
	// connecting vector from particle i to j
	std::vector<double> r_ij = {0.,0.,0.}; 
	double rSquared = 0.8; // distance of particles i,j squared
	double radiiSquared = 1.; // squared sum of particles i,j radii
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
			force->repulsionForce(
				forceI,
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
		}
	}
}

/* The reversible version of this function also calculates energies.*/
void Simulation::calculateRepulsionForcesRev()
{
	std::vector<double> forceI = {0.,0.,0.};
	std::vector<double> forceJ = {0.,0.,0.};
	// connecting vector from particle i to j
	std::vector<double> r_ij   = {0.,0.,0.}; 
	double rSquared     = 0.8; // distance of particles i,j squared
	double radiiSquared = 1.; // squared sum of particles i,j radii
	double energy       = 0.; // interaction energy of particle pair (i,j)
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
				energy,
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
			activeParticles[i].singleEnergy += energy;
			activeParticles[j].singleEnergy += energy;
		}
	}
}

void Simulation::updateOldPositions()
{
	for (int i=0; i<activeParticles.size(); i++) {
		activeParticles[i].oldPosition = activeParticles[i].position;
	}
}

void Simulation::resetSingleEnergies()
{
	for (int i=0; i<activeParticles.size(); i++) {
		activeParticles[i].singleEnergy = 0.;
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
