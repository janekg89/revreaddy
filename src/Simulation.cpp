/* Simulation.cpp
 * author: Christoph Froehner
 */

#include "Simulation.h"

Simulation::Simulation()//bool hasDefaultTypes)
{
	this->random            = new Random("ranlxs0");
	this->typeDict          = new TypeDict();
	this->timestep          = 0.001;
	this->cumulativeRuntime = 0.;
	this->temperature       = 1.;
	this->kBoltzmann        = 1.;
	this->isPeriodic        = true;
	this->boxsize           = 10.;
	this->energy            = 0.;
	this->oldEnergy         = 0.;
	this->currentAcceptance = 1.;
	this->acceptions        = 0;
	this->rejections        = 0;
	this->isReversible      = true;
	this->uniqueIdCounter   = 0;
	this->useNeighborList   = true;
}

Simulation::~Simulation()
{
	delete this->random;
	delete this->typeDict;
}

void Simulation::run()
{
	std::cout << "Simulation started at time: "
	          << this->cumulativeRuntime << "\n";
	this->resetForces();
	this->energy = 0.;
	this->calculateInteractionForcesEnergies();
	this->calculateGeometryForcesEnergies();
	this->recordObservables(0);
	for (unsigned long int timeIndex = 1; timeIndex < maxTime; timeIndex++)
	{
		// groupForces()
		this->saveOldState();
		this->propagate(); // propose
		this->resetForces();
		this->energy = 0.;
		if (this->useNeighborList) {
			this->calculateInteractionForcesEnergies(); // calculate energy and force
		}
		else {
			this->calculateInteractionForcesEnergiesNaive();
		}
		this->calculateGeometryForcesEnergies();
		if (this->isReversible) { this->acceptOrReject(); }
		this->cumulativeRuntime += this->timestep;
		this->recordObservables(timeIndex);
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

void Simulation::recordObservables(unsigned long int timeIndex)
{
	for (auto* obs : this->observables) {
		if (timeIndex % obs->recPeriod == 0) {
			obs->record(this->activeParticles, this->cumulativeRuntime);
		}
	}
}

void Simulation::calculateInteractionForcesEnergies()
{
	double minimalLength = this->boxsize;
	for (unsigned int i=0; i<possibleForces.size(); i++) {
		if (possibleForces[i]->cutoff < minimalLength) {
			minimalLength = possibleForces[i]->cutoff;
		}
	}
	double counter = 1.;
	unsigned int numberBoxes = 0;
	while ( (this->boxsize / counter) > minimalLength) {
		numberBoxes += 1;
		counter += 1;
	}
	/* if n=3 we will have 9 subboxes of length L/n, which
	 * will result in having to check every box. This is
	 * as inefficient as double looping. So:
	 * ONLY construct neighborlist if we have at least 16
	 * subboxes or n>3 */
	if ( numberBoxes > 3 ) {
		this->calculateInteractionForcesEnergiesWithLattice(numberBoxes);
	}
	else {
		this->calculateInteractionForcesEnergiesNaive();
	}
}

void Simulation::calculateInteractionForcesEnergiesNaive()
{
	for (unsigned int i=0; i<activeParticles.size(); i++) {
		for (unsigned int j=i+1; j<activeParticles.size(); j++) {
			this->calculateSingleForceEnergy(i, j);
		}
	}
}

void Simulation::calculateInteractionForcesEnergiesWithLattice(unsigned int numberBoxes)
{
	// construct neighborlist
	double n = (double) numberBoxes;
	double boxLength = this->boxsize / n;

	std::vector< std::vector< std::vector< std::vector<unsigned int> > > >
	neighborList(numberBoxes,
		std::vector< std::vector< std::vector<unsigned int> > >(numberBoxes,
			std::vector< std::vector<unsigned int> >(numberBoxes,
				std::vector<unsigned int>(0) ) ) );

	// reserve space for numberBoxes^3 lists containing particle indices
	neighborList.resize(numberBoxes);
	for (unsigned int x; x<numberBoxes; x++) {
		neighborList[x].resize(numberBoxes);
		for (unsigned int y; y<numberBoxes; y++) {
			neighborList[x][y].resize(numberBoxes);
		}
	}
	double delX = 0.;
	double delY = 0.;
	double delZ = 0.;
	unsigned int xIndex = 0;
	unsigned int yIndex = 0;
	unsigned int zIndex = 0;
	// find the right box triplet [x][y][z] for each particle
	for (unsigned int j=0; j<activeParticles.size(); j++) {
		delX = activeParticles[j].position[0] + 0.5*this->boxsize;
		xIndex = (unsigned int) floor(delX / boxLength);
		delY = activeParticles[j].position[1] + 0.5*this->boxsize;
		yIndex = (unsigned int) floor(delY / boxLength);
		delZ = activeParticles[j].position[2] + 0.5*this->boxsize;
		zIndex = (unsigned int) floor(delZ / boxLength);
		// add the particles index to the list of the corresponding box
		neighborList[xIndex][yIndex][zIndex].push_back(j);
	}
	// neighborList created
	// set up vector NxN filled with bools (false initially)
	std::vector< std::vector<bool> > alreadyCalculatedPairs;
	alreadyCalculatedPairs.resize(activeParticles.size());
	for (unsigned int i=0; i<alreadyCalculatedPairs.size(); i++) {
		alreadyCalculatedPairs[i].resize(activeParticles.size());
		std::fill(
			alreadyCalculatedPairs[i].begin(),
			alreadyCalculatedPairs[i].end(),
			false);
	}
	signed int otherX = 0;
	signed int otherY = 0;
	signed int otherZ = 0;
	for (unsigned int x=0; x<numberBoxes; x++)
	for (unsigned int y=0; y<numberBoxes; y++)
	for (unsigned int z=0; z<numberBoxes; z++) {
		for (signed int x_i = -1; x_i < 2; x_i++)
		for (signed int y_i = -1; y_i < 2; y_i++)
		for (signed int z_i = -1; z_i < 2; z_i++) {
			if ( (x_i==0) && (y_i==0) && (z_i==0) ) {
				for (unsigned int i=0;   i<neighborList[x][y][z].size(); i++)
				for (unsigned int j=i+1; j<neighborList[x][y][z].size(); j++) {
					this->calculateSingleForceEnergy(
						neighborList[x][y][z][i],
						neighborList[x][y][z][j]);
					alreadyCalculatedPairs
						[neighborList[x][y][z][i]]
						[neighborList[x][y][z][j]] = true;
					alreadyCalculatedPairs
						[neighborList[x][y][z][j]]
						[neighborList[x][y][z][i]] = true;
				}
			}
			else {
				//determine the "other" subbox. detect "over/underflows"
				otherX = x + x_i;
				if (otherX == -1) {otherX = numberBoxes - 1;}
				if (otherX == numberBoxes) {otherX = 0;}
				otherY = y + y_i;
				if (otherY == -1) {otherY = numberBoxes - 1;}
				if (otherY == numberBoxes) {otherY = 0;}
				otherZ = z + z_i;
				if (otherZ == -1) {otherZ = numberBoxes - 1;}
				if (otherZ == numberBoxes) {otherZ = 0;}
				for (unsigned int i=0; i<neighborList[x][y][z].size(); i++)
				for (unsigned int j=0; j<neighborList[otherX][otherY][otherZ].size(); j++) {
					// if not already calculated
					if (! alreadyCalculatedPairs
						[ neighborList[x][y][z][i] ]
						[ neighborList[otherX][otherY][otherZ][j] ] ) {
						// calculate interaction
						this->calculateSingleForceEnergy(
							neighborList[x][y][z][i],
							neighborList[otherX][otherY][otherZ][j] );
						alreadyCalculatedPairs
							[ neighborList[x][y][z][i] ]
							[ neighborList[otherX][otherY][otherZ][j] ] = true;
						alreadyCalculatedPairs
							[ neighborList[otherX][otherY][otherZ][j] ]
							[ neighborList[x][y][z][i] ] = true;
					}
				}
			}
		}
	}
}

void Simulation::calculateSingleForceEnergy(
	unsigned int indexI,
	unsigned int indexJ)
{
	std::vector<double> forceI = {0.,0.,0.};
	std::vector<double> forceJ = {0.,0.,0.};
	// interaction energy of particle pair (i,j)
	double energyBuffer = 0.; 
	// look for force that affects the pair (i,j)
	for (unsigned int k=0; k<this->possibleForces.size(); k++) {
		if (
			this->possibleForces[k]->isAffected(
				this->activeParticles[indexI].typeId,
				this->activeParticles[indexJ].typeId)
		) {
			// connecting vector from particle i to j
			std::vector<double> r_ij = {0.,0.,0.}; 
			getMinDistanceVector(
				r_ij,
				activeParticles[indexI].position, 
				activeParticles[indexJ].position, 
				this->isPeriodic, 
				this->boxsize);
			// distance of particle i,j squared
			double rSquared = r_ij[0]*r_ij[0] + r_ij[1]*r_ij[1] + r_ij[2]*r_ij[2];
			// radius of particle i
			double radiusI = this->typeDict->radii[activeParticles[indexI].typeId]; 
			// radius of particle j
			double radiusJ = this->typeDict->radii[activeParticles[indexJ].typeId];
			// squared sum of particles i,j radii
			double radiiSquared = pow(radiusI + radiusJ, 2.);
			// actual force call here
			this->possibleForces[k]->calculateForceEnergy(
				forceI,//out
				energyBuffer,
				r_ij,//in
				rSquared,
				radiiSquared);
			forceJ[0] = -1. * forceI[0];
			forceJ[1] = -1. * forceI[1];
			forceJ[2] = -1. * forceI[2];
			activeParticles[indexI].addForce(forceI);
			activeParticles[indexJ].addForce(forceJ);
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
	this->currentAcceptance = acceptance;
	// accept or reject
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

// TODO check if all numbers are > 0.
void Simulation::new_Type(
	std::string name,
	double radius,
	double diffusionConstant,
	double reactionRadius)
{
	this->typeDict->newType(
		name,
		radius, 
		diffusionConstant,
		reactionRadius);
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

void Simulation::writeLastObservableToFile()
{
	if (this->observables.size() > 0) {
		this->observables.back()->writeBufferToFile();
	}
	else {
		std::cout << "Error: There are no observables to write" << std::endl;
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

void Simulation::deleteLastObservable()
{
	/* first delete the observable, then its pointer in the vector */
	if (this->observables.size() > 0) {
		delete this->observables.back();
		this->observables.pop_back();
	}
}

void Simulation::new_Trajectory(unsigned long int recPeriod, std::string filename)
{
	Trajectory * obs = new Trajectory(filename);
	obs->recPeriod = recPeriod;
	this->observables.push_back(obs);
}

void Simulation::new_RadialDistribution(
	unsigned long int recPeriod,
	std::string filename,
	std::vector<double> ranges,
	std::vector< std::vector<unsigned int> > considered)
{
	RadialDistribution * rad = new RadialDistribution(
		ranges,
		this->isPeriodic,
		this->boxsize,
		considered,
		filename);
	rad->recPeriod = recPeriod;
	this->observables.push_back(rad);
}

void Simulation::new_MeanSquaredDisplacement(
	unsigned long int recPeriod,
	std::string filename,
	unsigned int particleTypeId)
{
	MeanSquaredDisplacement * msd = new MeanSquaredDisplacement(
		this->activeParticles,
		particleTypeId,
		this->cumulativeRuntime,
		this->boxsize,
		filename);
	msd->recPeriod = recPeriod;
	this->observables.push_back(msd);
}

void Simulation::new_ProbabilityDensity(
	unsigned long int recPeriod,
	std::string filename,
	unsigned int pTypeId,
	std::vector<double> range,
	unsigned int coord)
{
	ProbabilityDensity * prob = new ProbabilityDensity(
		this->activeParticles,
		pTypeId,
		range,
		coord,
		filename);
	prob->recPeriod = recPeriod;
	this->observables.push_back(prob);
}

void Simulation::new_Energy(unsigned long int recPeriod, std::string filename)
{
	Energy * ener = new Energy(
		recPeriod,
		0,
		this,
		filename);
	this->observables.push_back(ener);
}

void Simulation::new_Acceptance(unsigned long int recPeriod, std::string filename)
{
	Acceptance * acc = new Acceptance(
		recPeriod,
		0,
		this,
		filename);
	this->observables.push_back(acc);
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

void Simulation::deleteAllForces()
{
	/* first delete the forces, since they we're allocated with 'new'
	 * then erase the pointers from the vector */
	for (auto* f : this->possibleForces) {
		delete f;
	}
	this->possibleForces.erase(
		this->possibleForces.begin(),
		this->possibleForces.begin() + this->possibleForces.size()
	);
}

void Simulation::new_SoftRepulsion(
	std::string name,
	std::vector<unsigned int> affectedTuple,
	double repulsionStrength)
{
	if (affectedTuple.size() != 2) {
		std::cout << "Error: The given tuple must be of length 2\n";
		return;
	}
	if ( (affectedTuple[0] > ( this->typeDict->names.size() - 1) ) 
	  || (affectedTuple[1] > ( this->typeDict->names.size() - 1) ) ) {
		std::cout << "Error: The given particle type(s) do not exist. "
		          << "Make sure to add them first\n";
		return;
	}
	if ( repulsionStrength <= 0. ) {
		std::cout << "Error: The repulsion strength must be larger than zero\n";
		return;
	}
	SoftRepulsion * soft = new SoftRepulsion(
		name,
		affectedTuple,
		repulsionStrength);
	// set cutoff correctly
	soft->cutoff = this->typeDict->radii[affectedTuple[0]] 
	             + this->typeDict->radii[affectedTuple[1]];
	this->possibleForces.push_back(soft);
	std::cout << "Info: SoftRepulsion interaction added to possibleForces\n";
}

void Simulation::new_LennardJones(
	std::string name,
	std::vector<unsigned int> affectedTuple,
	double epsilon)
{
	if (affectedTuple.size() != 2) {
		std::cout << "Error: The given tuple must be of length 2" << std::endl;
		return;
	}
	if ( (affectedTuple[0] > ( this->typeDict->names.size() - 1) ) 
	  || (affectedTuple[1] > ( this->typeDict->names.size() - 1) ) ) {
		std::cout << "Error: The given particle type(s) do not exist. "
		          << "Make sure to add them first\n";
		return;
	}
	if ( epsilon <= 0. ) {
		std::cout << "Error: The given epsilon must be larger than zero\n";
		return;
	}
	LennardJones * lj = new LennardJones(
		name,
		affectedTuple,
		epsilon);
	// set cutoff correctly
	lj->cutoff = 2.5 * ( this->typeDict->radii[affectedTuple[0]]
	                   + this->typeDict->radii[affectedTuple[1]] );
	this->possibleForces.push_back(lj);
	std::cout << "Info: LennardJones interaction added to possibleForces\n";
}

unsigned int Simulation::getNumberForces()
{
	return this->possibleForces.size();
}

std::string Simulation::getForceName(unsigned int i)
{
	return this->possibleForces[i]->name;
}

std::string Simulation::getForceType(unsigned int i)
{
	return this->possibleForces[i]->type;
}

std::vector<unsigned int> Simulation::getForceAffectedTuple(unsigned int i)
{
	return this->possibleForces[i]->affectedTuple;
}

std::vector<double> Simulation::getForceParameters(unsigned int i)
{
	return this->possibleForces[i]->parameters;
}
