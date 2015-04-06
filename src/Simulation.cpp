/* Simulation.cpp
 * author: Christoph Froehner
 */

#include "Simulation.h"

Simulation::Simulation(bool hasDefaultTypes)
{
	this->random            = new Random("ranlxs0");
	this->force             = new Force();
	this->typeDict          = new TypeDict();
	this->timestep          = 0.001;
	this->cumulativeRuntime = 0.;
	this->temperature       = 1.;
	this->kBoltzmann        = 1.;
	this->repulsionStrength = 1.;
	this->isPeriodic        = true;
	this->boxsize           = 10.;
	this->energy            = 0.;
	this->oldEnergy         = 0.;
	this->acceptions        = 0;
	this->rejections        = 0;
	this->isReversible      = true;
	this->uniqueIdCounter   = 0;
	if (hasDefaultTypes) {
		this->typeDict->newType("default", 1., 0., 1., 0); // 0
		this->typeDict->newType("lj", 1., 1., 1., 1); // 1
		this->typeDict->newType("soft", 1., 1., 1., 2); // 2
	}
}

Simulation::~Simulation()
{
	delete random;
	delete force;
}

void Simulation::run()
{
	std::cout << "Simulation started at time: "
	          << this->cumulativeRuntime << "\n";
	this->resetForces();
	this->energy = 0.;
	this->calculateRepulsionForcesEnergies();
	this->calculateGeometryForcesEnergies();
	this->recordObservables(this->cumulativeRuntime);
	for (unsigned long int t = 1; t < maxTime; t++)
	{
		// groupForces()
		this->saveOldState();
		this->propagate(); // propose
		this->resetForces();
		this->energy = 0.;
		this->calculateRepulsionForcesEnergies(); // calculate energy and force
		this->calculateGeometryForcesEnergies();
		if (this->isReversible) { this->acceptOrReject(); }
		this->cumulativeRuntime += this->timestep;
		this->recordObservables(this->cumulativeRuntime);
	}
	std::cout << "Simulation finished at time: "
	          << this->cumulativeRuntime << "\n";
}

void Simulation::saveOldState()
{
	this->oldEnergy = this->energy;
	for (int i=0; i<activeParticles.size(); i++) {
		activeParticles[i].oldPosition 
			= activeParticles[i].position;
		activeParticles[i].oldBoxCoordinates 
			= activeParticles[i].boxCoordinates;
		activeParticles[i].oldForce 
			= activeParticles[i].cumulativeForce;
	}
}

void Simulation::propagate()
{
	std::vector<double> noiseTerm = {0.,0.,0.};
	std::vector<double> forceTerm = {0.,0.,0.};
	double noisePrefactor = 1.;
	double forcePrefactor = 1.;
	double diffConst = 1.; //diffusion constant of current particle
	for (int i=0; i<activeParticles.size(); i++)
	{
		// look up particles' diffusion constant from its typeId
		diffConst=this->typeDict->diffusionConstants[activeParticles[i].typeId];

		noiseTerm = random->normal3D();
		noisePrefactor = sqrt(2. * diffConst * timestep);
		noiseTerm[0] *= noisePrefactor;
		noiseTerm[1] *= noisePrefactor;
		noiseTerm[2] *= noisePrefactor;

		forcePrefactor = timestep * diffConst / (kBoltzmann * temperature);
		forceTerm[0] = activeParticles[i].cumulativeForce[0] * forcePrefactor;
		forceTerm[1] = activeParticles[i].cumulativeForce[1] * forcePrefactor;
		forceTerm[2] = activeParticles[i].cumulativeForce[2] * forcePrefactor;

		activeParticles[i].move(noiseTerm);
		activeParticles[i].move(forceTerm);

		if (isPeriodic)
		{
			if (activeParticles[i].position[0] < (-0.5 * boxsize) ) {
				activeParticles[i].position[0] += boxsize;
				activeParticles[i].boxCoordinates[0] -= 1;
			}
			else if (activeParticles[i].position[0] >= (0.5 * boxsize) ) {
				activeParticles[i].position[0] -= boxsize;
				activeParticles[i].boxCoordinates[0] += 1;
			}
			if (activeParticles[i].position[1] < (-0.5 * boxsize) ) {
				activeParticles[i].position[1] += boxsize;
				activeParticles[i].boxCoordinates[1] -= 1;
			}
			else if (activeParticles[i].position[1] >= (0.5 * boxsize) ) {
				activeParticles[i].position[1] -= boxsize;
				activeParticles[i].boxCoordinates[1] += 1;
			}
			if (activeParticles[i].position[2] < (-0.5 * boxsize) ) {
				activeParticles[i].position[2] += boxsize;
				activeParticles[i].boxCoordinates[2] -= 1;
			}
			else if (activeParticles[i].position[2] >= (0.5 * boxsize) ) {
				activeParticles[i].position[2] -= boxsize;
				activeParticles[i].boxCoordinates[2] += 1;
			}
		}
	}
}

void Simulation::recordObservables(double t)
{
	for (auto* obs : this->observables) {
		obs->record(this->activeParticles, t);
	}
}

void Simulation::calculateRepulsionForcesEnergies()
{
	std::vector<double> forceI = {0.,0.,0.};
	std::vector<double> forceJ = {0.,0.,0.};
	// connecting vector from particle i to j
	std::vector<double> r_ij   = {0.,0.,0.}; 
	double rSquared     = 0.8; // distance of particles i,j squared
	double radiiSquared = 1.; // squared sum of particles i,j radii
	double energyBuffer = 0.; // interaction energy of particle pair (i,j)
	double radiusI = 0.; // radius of particle i
	double radiusJ = 0.; // radius of particle j
	for (int i=0; i<activeParticles.size(); i++) {
		radiusI = this->typeDict->radii[activeParticles[i].typeId];
		for (int j=i+1; j<activeParticles.size(); j++) {
			radiusJ = this->typeDict->radii[activeParticles[j].typeId];
			getMinDistanceVector(
				r_ij,
				activeParticles[i].position, 
				activeParticles[j].position, 
				this->isPeriodic, 
				this->boxsize);
			rSquared = r_ij[0]*r_ij[0] + r_ij[1]*r_ij[1] + r_ij[2]*r_ij[2];
			radiiSquared = pow(radiusI + radiusJ, 2.);
			force->repulsionForceEnergy(
				forceI,
				energyBuffer,
				r_ij,
				rSquared,
				radiiSquared,
				repulsionStrength, 
				this->typeDict->forceTypes[activeParticles[i].typeId],
				this->typeDict->forceTypes[activeParticles[j].typeId]);
			forceJ[0] = -1. * forceI[0];
			forceJ[1] = -1. * forceI[1];
			forceJ[2] = -1. * forceI[2];
			activeParticles[i].addForce(forceI);
			activeParticles[j].addForce(forceJ);
			this->energy += energyBuffer;
		}
	}
}

void Simulation::calculateGeometryForcesEnergies()
{
	std::vector<double> forceI = {0.,0.,0.};
	double energyBuffer = 0.;
	for (int i=0; i<activeParticles.size(); i++) {
		for (int j=0; j<geometries.size(); j++) {
			if (geometries[j]->doesInteract(activeParticles[i].typeId)) {
				geometries[j]->forceEnergy(
					forceI,
					energyBuffer,
					activeParticles[i].position,
					this->typeDict->radii[activeParticles[i].typeId]);
				activeParticles[i].addForce(forceI);
				this->energy += energyBuffer;
			}
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
		secondTerm += 
		          this->typeDict->diffusionConstants[activeParticles[i].typeId]
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
		else {/* reject = restore old positions,boxCoordinates,forces,energy */
			rejections += 1;
			this->energy = this->oldEnergy;
			for (int i=0; i<activeParticles.size(); i++) {
				activeParticles[i].position 
					= activeParticles[i].oldPosition;
				activeParticles[i].cumulativeForce
					= activeParticles[i].oldForce;
				activeParticles[i].boxCoordinates
					= activeParticles[i].oldBoxCoordinates;
			}
		}
	}
}

void Simulation::addParticle(
	std::vector<double> initPos,
	unsigned int particleTypeId)
{
	if (particleTypeId >= this->typeDict->getNumberOfTypes()) {
		std::cout << "The given particle type does not exist!\n"
		          << "Particle is not created.\n";
		return;
	}
	Particle * particle = new Particle();
	if ( initPos.size() == 3 ) { particle->position  = initPos; }
	else {
		std::cout << "Particles' initial position has dimension mismatch!\n" 
		          << "Particle will be placed at {0,0,0}.\n";	
		particle->position = {0., 0., 0.};
	}
	particle->typeId = particleTypeId;
	particle->uniqueId = this->uniqueIdCounter;
	this->uniqueIdCounter += 1;
	this->activeParticles.push_back(*particle);//push_back copies arg into vec
	delete particle;
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
		std::cout << "New position has dimension mismatch!\n"
		          << "Particle remains at its old position.\n";
	}
}

unsigned int Simulation::getTypeId(int index) {
	return this->activeParticles[index].typeId;
}

void Simulation::setTypeId(int index, unsigned int typeId) 
{
	if (typeId >= this->typeDict->getNumberOfTypes()) {
		std::cout << "The given particle type does not exist!\n"
		          << "Particle is not created.\n";
		return;
	}
	this->activeParticles[index].typeId = typeId;
}

void Simulation::new_Type(
	std::string name,
	double radius,
	double diffusionConstant,
	double reactionRadius,
	unsigned int forceType) 
{
	this->typeDict->newType(
		name,
		radius, 
		diffusionConstant,
		reactionRadius,
		forceType);
}

std::vector<std::string> Simulation::getDictNames() {
	return this->typeDict->names;
}

std::vector<double> Simulation::getDictRadii() {
	return this->typeDict->radii;
}

std::vector<double> Simulation::getDictDiffusionConstants() {
	return this->typeDict->diffusionConstants;
}

std::vector<double> Simulation::getDictReactionRadii() {
	return this->typeDict->reactionRadii;
}

std::vector<unsigned int> Simulation::getDictForceTypes() {
	return this->typeDict->forceTypes;
}

int Simulation::getParticleNumber() {
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
	this->observables.push_back(obs);
}

void Simulation::new_RadialDistribution(
	std::string filename,
	std::vector<double> ranges,
	std::vector< std::vector<unsigned int> > considered)
{
	RadialDistribution * rad = new RadialDistribution(
		ranges,
		this->isPeriodic,
		this->boxsize,
		considered);
	rad->filename = filename;
	this->observables.push_back(rad);
}

void Simulation::new_MeanSquaredDisplacement(
	std::string filename,
	unsigned int particleTypeId)
{
	MeanSquaredDisplacement * msd = new MeanSquaredDisplacement(
		this->activeParticles,
		particleTypeId,
		this->cumulativeRuntime,
		this->boxsize);
	msd->filename = filename;
	this->observables.push_back(msd);
}

void Simulation::new_ProbabilityDensity(
	std::string filename,
	unsigned int pTypeId,
	std::vector<double> range,
	unsigned int coord)
{
	ProbabilityDensity * prob = new ProbabilityDensity(
		this->activeParticles,
		pTypeId,
		range,
		coord);
	prob->filename = filename;
	this->observables.push_back(prob);
}

void Simulation::deleteAllGeometries()
{
	/* first delete the geometries, since they we're allocated with 'new'
	 * then erase the pointers from the vector */
	for (auto* geo : this->geometries) {
		delete geo;
	}
	this->geometries.erase(
		this->geometries.begin(),
		this->geometries.begin() + this->geometries.size()
	);
}

void Simulation::new_Wall(
	std::vector<double> normal,
	std::vector<double> point,
	double strength,
	std::vector<unsigned int>& particleTypeIds)
{
	Wall * wall = new Wall(
		normal,
		point,
		strength,
		particleTypeIds);
	this->geometries.push_back(wall);
}

void Simulation::new_DoubleWellZ(
	double distanceMinima,
	double strength,
	std::vector<unsigned int> particleTypeIds)
{
	DoubleWellZ * well = new DoubleWellZ(
		distanceMinima,
		strength,
		particleTypeIds);
	this->geometries.push_back(well);
}
