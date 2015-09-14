/* Simulation.cpp
 * author: Christoph Froehner
 */

#include "Simulation.h"

Simulation::Simulation()
{
	this->random                = new Random("ranlxs0");
	//this->typeDict              = new TypeDict();
	this->timestep              = 0.001;
	this->cumulativeRuntime     = 0.;
	this->temperature           = 1.;
	this->kBoltzmann            = 1.;
	this->isPeriodic            = true;
	this->boxsize               = 10.;
	this->energy                = 0.;
	this->oldEnergy             = 0.;
	this->acceptProbDynamics    = 1.;
	this->acceptProbReactions   = 1.;
	this->acceptionsDynamics    = 0;
	this->rejectionsDynamics    = 0;
	this->isReversibleDynamics  = true;
	this->isReversibleReactions = true;
	this->uniqueIdCounter       = 0;
	this->useNeighborList       = true;
	this->reactionPropagation   = 0;
}

Simulation::~Simulation()
{
	delete this->random;
}

/* TODO after one run() the observables' files should 
 * be closed so that the intermediate results can be 
 * accessed (e.g. from python). is that so??? 
 * When run is called the observables' files should 
 * be opened again.*/
void Simulation::run()
{
	std::cout << "Info: Started at simulation-time: "
	          << this->cumulativeRuntime << std::endl;
	std::cout << "Info: Run ..." << std::endl;
	clock_t timer = clock();
	this->resetForces();
	this->resetActivePairs();
	this->energy = 0.;
	this->calculateInteractionForcesEnergies();
	this->calculateGeometryForcesEnergies();
	this->recordObservables(0);
	double acceptance = 1.;
	bool isStepAccepted = true;
	for (unsigned long int timeIndex = 1; timeIndex < maxTime; timeIndex++)
	{
		/* Dynamics */
		this->saveOldState();
		this->propagateDynamics(); // propose
		this->resetForces();
		this->resetActivePairs();
		this->energy = 0.;
		this->calculateInteractionForcesEnergies(); // calculate energy and force
		this->calculateGeometryForcesEnergies();
		acceptance = this->acceptanceDynamics();
		this->acceptProbDynamics = acceptance;
		isStepAccepted = this->acceptOrReject(acceptance);
		if (this->isReversibleDynamics && ( ! isStepAccepted ) ) {
			this->restoreOldState();
			this->rejectionsDynamics += 1;
		}
		else { this->acceptionsDynamics += 1; }

		/* Reactions */
		this->saveOldState();
		acceptance = this->propagateReactions();
		this->resetForces();
		this->resetActivePairs();
		this->energy = 0.;
		this->calculateInteractionForcesEnergies();
		this->calculateGeometryForcesEnergies();
		acceptance *= this->acceptanceReactions();
		this->acceptProbReactions = acceptance;
		isStepAccepted = this->acceptOrReject(acceptance);
		if (this->isReversibleReactions && ( ! isStepAccepted ) ) {
			this->restoreOldState();
			this->rejectionsReactions += 1;
		}
		else { this->acceptionsReactions += 1; }

		/* Advance clock */
		this->cumulativeRuntime += this->timestep;
		this->recordObservables(timeIndex);
	}
	timer = clock() - timer;
	std::cout << "Info: Finished at simulation-time: "
	          << this->cumulativeRuntime << std::endl;
	std::cout << "Info: Needed computation-time for run(): "
	          << ((float) timer) / CLOCKS_PER_SEC 
	          << " s" << std::endl;
}

void Simulation::saveOldState()
{
	this->oldEnergy          = this->energy;
	this->oldActiveParticles = this->activeParticles;
	this->oldActivePairs     = this->activePairs;
}

void Simulation::restoreOldState()
{
	this->energy          = this->oldEnergy;
	this->activeParticles = this->oldActiveParticles;
	this->activePairs     = this->oldActivePairs;
}

void Simulation::propagateDynamics()
{
	std::vector<double> noiseTerm = {0.,0.,0.};
	std::vector<double> forceTerm = {0.,0.,0.};
	double noisePrefactor = 1.;
	double forcePrefactor = 1.;
	double diffConst = 1.; //diffusion constant of current particle
	for (unsigned long int i=0; i<activeParticles.size(); i++)
	{
		// look up particles' diffusion constant from its typeId
		diffConst=this->typeDict[activeParticles[i].typeId].diffusionConstant;

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

//TODO
double Simulation::propagateReactions()
{
	/* First construct the queue reactionCandidates, then shuffle it
	 * and perform the reactions subsequently. If one particle
	 * cannot be found via its uniqueId, because it was destroyed in
	 * a previous reaction, skip this event and check the next 
	 * candidate. Note that a ReactionEvent will be performed with
	 * its predefined probability. */
	std::vector<ReactionEvent> reactionCandidates;
	/* Find bimolecular candidates */
	for (unsigned long i=0; i<this->activePairs.size(); i++) {
		// extract indices of the pair
		unsigned long particle1Index = activePairs[i][0];
		unsigned long particle2Index = activePairs[i][1];
		// determine types of the pair
		unsigned int particle1Type
			= activeParticles[particle1Index].typeId;
		unsigned int particle2Type
			= activeParticles[particle2Index].typeId;
		std::vector<unsigned int> types = {particle1Type, particle2Type};
		for (unsigned int j=0; j<this->possibleReactions; j++) {
			if ( possibleReactions[j]->isAffectedForward(types) ) {
				std::vector<unsigned long long> participants;
				participants.push_back(
					activeParticles[particle1Index].uniqueId);
				participants.push_back(
					activeParticles[particle2Index].uniqueId);
				ReactionEvent event(
					j, // reactionId
					true, // forwardOrBackward
					participants); // uniqueIds of reaction participants
			}
			else if ( possibleReaction[j]->isAffectedBackward(types) ) {
				std::vector<unsigned long long> participants;
				participants.push_back(
					activeParticles[particle1Index].uniqueId);
				participants.push_back(
					activeParticles[particle2Index].uniqueId);
				ReactionEvent event(
					j, // reactionId
					false, // forwardOrBackward
					participants); // uniqueIds of reaction participants			
			}
		}
	}
	/* Find unimolecular candidates */
	for (unsigned long i=0; i<this->activeParticles.size(), i++) {
		// TODO
	}
	std::random_shuffle(reactionCandidates.begin(), reactionCandidates.end());
	/* Perform reactions subsequently */
	for (unsigned long i=0; i<reactionCandidates.size(); i++) {
		// TODO
	}
	double acceptance = 1.;
	return acceptance;
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
	/* We want to find out if neighborlist pays off. 
	 * First determine the minimal size of subboxes,
	 * therefore check all interaction distances and
	 * all reaction radii. */
	double minimalLength = this->boxsize;
	/* check interaction distances */
	for (unsigned int i=0; i<possibleInteractions.size(); i++) {
		if (possibleInteractions[i]->cutoff < minimalLength) {
			minimalLength = possibleInteractions[i]->cutoff;
		}
	}
	/* check reaction radii combinations (i,j) with i <= j */
	for (unsigned int i=0; i<this->typeDict.size(); i++) {
		for (unsigned int j=i; j<this->typeDict.size(); j++) {
			double reactionRadii
				= typeDict[i].reactionRadius + typeDict[j].reactionRadius;
			if (reactionRadii < minimalLength) {
				minimalLength = reactionRadii;
			}
		}
	}
	double counter = 1.;
	unsigned int numberBoxes = 0;
	while ( (this->boxsize / counter) > minimalLength) {
		numberBoxes += 1;
		counter += 1.;
	}
	/* if n=3 we will have 9 subboxes of length L/n, which
	 * will result in having to check every box. This is
	 * as inefficient as double looping. So:
	 * ONLY construct neighborlist if we have at least 16
	 * subboxes or n>3 */
	if ( ( numberBoxes > 3 ) && this->useNeighborList ) {
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

void Simulation::calculateInteractionForcesEnergiesWithLattice(
	unsigned int numberBoxes)
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
	/* First gather all parameters, distances and reaction
	 * radii of particles, check if they are within reaction
	 * distance and then check for forces that apply to the given
	 * pair of particles. */
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
	// radii of particle i and j
	double radiusI = this->typeDict[activeParticles[indexI].typeId].radius; 
	double radiusJ = this->typeDict[activeParticles[indexJ].typeId].radius;
	// squared sum of particles i,j radii
	double radiiSquared = pow(radiusI + radiusJ, 2.);
	// squared sum of reactionRadii of particle i and j
	double reactionRadiiSquared
		= pow(
		  this->typeDict[activeParticles[indexI].typeId].reactionRadius
		+ this->typeDict[activeParticles[indexJ].typeId].reactionRadius
		, 2.);
	// check if particles i, j are within reactive distance
	// if so, their unique ids will be added to activePairs
	if ( rSquared <= reactionRadiiSquared ) {
		std::vector<unsigned long> activePair;
		activePair.push_back(indexI);
		activePair.push_back(indexJ);
		activePair.shrink_to_fit();
		this->activePairs.push_back(activePair);
	}
	// look for force that affects the pair (i,j)
	for (unsigned int k=0; k<this->possibleInteractions.size(); k++) {
		if (
			this->possibleInteractions[k]->isAffected(
				this->activeParticles[indexI].typeId,
				this->activeParticles[indexJ].typeId)
		) {
			// actual force call here
			this->possibleInteractions[k]->calculateForceEnergy(
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
					this->typeDict[activeParticles[i].typeId].radius);
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

void Simulation::resetActivePairs()
{
	this->activePairs.clear();
}

/* For the dynamical acceptance probability both states, old and new
 * have to have the same number of particles.  */
double Simulation::acceptanceDynamics()
{
	double acceptance = 1.;
	double firstTerm  = 0.;
	double secondTerm = 0.;
	std::vector<double> deltaX = {0.,0.,0.};
	for (int i=0; i<activeParticles.size(); i++) {
		getMinDistanceVector(
			deltaX,
			oldActiveParticles[i].position,
			activeParticles[i].position,
			this->isPeriodic,
			this->boxsize);
		firstTerm  += deltaX[0]
		            * ( oldActiveParticles[i].cumulativeForce[0]
		              + activeParticles[i].cumulativeForce[0] );
		firstTerm  += deltaX[1]
		            * ( oldActiveParticles[i].cumulativeForce[1]
		              + activeParticles[i].cumulativeForce[1] );
		firstTerm  += deltaX[2]
		            * ( oldActiveParticles[i].cumulativeForce[2]
		              + activeParticles[i].cumulativeForce[2] );
		secondTerm += 
		          this->typeDict[activeParticles[i].typeId].diffusionConstant
		            * ( activeParticles[i].cumulativeForce[0]
		              * activeParticles[i].cumulativeForce[0]
		              + activeParticles[i].cumulativeForce[1]
		              * activeParticles[i].cumulativeForce[1]
		              + activeParticles[i].cumulativeForce[2]
		              * activeParticles[i].cumulativeForce[2]
		              + oldActiveParticles[i].cumulativeForce[0]
		              * oldActiveParticles[i].cumulativeForce[0]
		              + oldActiveParticles[i].cumulativeForce[1]
		              * oldActiveParticles[i].cumulativeForce[1]
		              + oldActiveParticles[i].cumulativeForce[2]
		              * oldActiveParticles[i].cumulativeForce[2] );
	}
	firstTerm  *= 0.5;
	secondTerm *= this->timestep / (4. * this->kBoltzmann * this->temperature);
	acceptance = firstTerm + secondTerm + this->energy - this->oldEnergy;
	acceptance /= -1. * this->kBoltzmann * this->temperature;
	acceptance = exp( acceptance );
	return acceptance;
}

//TODO
double Simulation::acceptanceReactions()
{
	return 1.;
}

bool Simulation::acceptOrReject(double acceptance)
{
	if ( acceptance > 1. ) {
		/* accept */
		return true;
	}
	else {
		double uniform = this->random->uniform();
		if ( uniform < acceptance ) {
			/* accept */
			return true;
		}
		else {
			/* reject */
			return false;
		}
	}
}

long int Simulation::findParticleIndex(unsigned long long id)
{
	unsigned long int max = this->activeParticles.size() - 1;
	unsigned long int min = 0;
	unsigned long int mid = 0;
	while (max >= min) {
		mid = min + (max - min) / 2;
		if (this->activeParticles[mid].uniqueId == id) {return mid;}
		else if (this->activeParticles[mid].uniqueId < id) {min = mid + 1;}
		else {max = mid - 1;}
	}
	// particle was not found
	return -1;
}

void Simulation::addParticle(
	std::vector<double> initPos,
	unsigned int particleTypeId)
{
	if (particleTypeId >= this->typeDict.size() ) {
		std::cout << "Error: The given particle type does not exist!\n"
		          << "Particle is not created" << std::endl;
		return;
	}
	Particle particle;
	if ( initPos.size() == 3 ) { particle.position  = initPos; }
	else {
		std::cout << "Error: Particles' initial position has dimension mismatch!\n" 
		          << "Particle will be placed at {0,0,0}" << std::endl;	
		particle.position = {0., 0., 0.};
	}
	particle.typeId = particleTypeId;
	particle.uniqueId = this->uniqueIdCounter;
	this->uniqueIdCounter += 1;
	this->activeParticles.push_back(particle);//push_back copies arg into vec
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
		std::cout << "Error: New position has dimension mismatch!\n"
		          << "Particle remains at its old position" << std::endl;
	}
}

unsigned int Simulation::getTypeId(int index) {
	return this->activeParticles[index].typeId;
}

void Simulation::setTypeId(int index, unsigned int typeId) 
{
	if (typeId >= this->typeDict.size() ) {
		std::cout << "Error: The given particle type does not exist!\n"
		          << "Particle is not created" << std::endl;
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
	if ( radius < 0. ) {
		std::cout << "Error: The particle radius must be non-negative!\n"
		          << "New type is not created" << std::endl;
		return;
	}
	if ( diffusionConstant < 0. ) {
		std::cout << "Error: The diffusionConstant must be non-negative!\n"
		          << "New Type is not created" << std::endl;
		return;
	}
	if ( reactionRadius < 0. ) {
		std::cout << "Error: The reactionRadius must be non-negative!\n"
		          << "New Type is not created" << std::endl;
		return;
	}
	ParticleType pType(
		name,
		radius,
		diffusionConstant,
		reactionRadius);
	this->typeDict.push_back(pType);
}

unsigned int Simulation::getNumberOfTypes() {
	return this->typeDict.size();
}

std::string Simulation::getDictName(unsigned int i) {
	return this->typeDict[i].name;
}

double Simulation::getDictRadius(unsigned int i) {
	return this->typeDict[i].radius;
}

double Simulation::getDictDiffusionConstant(unsigned int i) {
	return this->typeDict[i].diffusionConstant;
}

double Simulation::getDictReactionRadius(unsigned int i) {
	return this->typeDict[i].reactionRadius;
}

unsigned int Simulation::getParticleNumber() {
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
	for (auto* f : this->possibleInteractions) {
		delete f;
	}
	this->possibleInteractions.erase(
		this->possibleInteractions.begin(),
		this->possibleInteractions.begin() + this->possibleInteractions.size()
	);
}

void Simulation::new_SoftRepulsion(
	std::string name,
	std::vector<unsigned int> affectedTuple,
	double repulsionStrength)
{
	if (affectedTuple.size() != 2) {
		std::cout << "Error: The given tuple must be of length 2" << std::endl;
		return;
	}
	if ( (affectedTuple[0] > ( this->typeDict.size() - 1) ) 
	  || (affectedTuple[1] > ( this->typeDict.size() - 1) ) ) {
		std::cout << "Error: The given particle type(s) do not exist.\n"
		          << "Make sure to add them first" << std::endl;
		return;
	}
	if ( repulsionStrength <= 0. ) {
		std::cout << "Error: The repulsion strength must be larger than zero"
		          << std::endl;
		return;
	}
	SoftRepulsion * soft = new SoftRepulsion(
		name,
		affectedTuple,
		repulsionStrength);
	// set cutoff correctly
	soft->cutoff = this->typeDict[affectedTuple[0]].radius 
	             + this->typeDict[affectedTuple[1]].radius;
	this->possibleInteractions.push_back(soft);
	std::cout << "Info: SoftRepulsion interaction added to possibleInteractions"
	          << std::endl;
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
	if ( (affectedTuple[0] > ( this->typeDict.size() - 1) ) 
	  || (affectedTuple[1] > ( this->typeDict.size() - 1) ) ) {
		std::cout << "Error: The given particle type(s) do not exist.\n"
		          << "Make sure to add them first" << std::endl;
		return;
	}
	if ( epsilon <= 0. ) {
		std::cout << "Error: The given epsilon must be larger than zero"
		          << std::endl;
		return;
	}
	LennardJones * lj = new LennardJones(
		name,
		affectedTuple,
		epsilon);
	// set cutoff correctly
	lj->cutoff = 2.5 * ( this->typeDict[affectedTuple[0]].radius
	                   + this->typeDict[affectedTuple[1]].radius );
	this->possibleInteractions.push_back(lj);
	std::cout << "Info: LennardJones interaction added to possibleInteractions"
	          << std::endl;
}

unsigned int Simulation::getNumberForces()
{
	return this->possibleInteractions.size();
}

std::string Simulation::getForceName(unsigned int i)
{
	return this->possibleInteractions[i]->name;
}

std::string Simulation::getForceType(unsigned int i)
{
	return this->possibleInteractions[i]->type;
}

std::vector<unsigned int> Simulation::getForceAffectedTuple(unsigned int i)
{
	return this->possibleInteractions[i]->affectedTuple;
}

std::vector<double> Simulation::getForceParameters(unsigned int i)
{
	return this->possibleInteractions[i]->parameters;
}
