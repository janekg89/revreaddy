/* Simulation.cpp
 * author: Christoph Froehner */
#define print(a) std::cout << a << std::endl;

#include "Simulation.h"

Simulation::Simulation()
{
	this->random = new Random("ranlxs0");
	this->world  = new World();
	this->config = new Config(this->world, this->random);
}

Simulation::~Simulation()
{
	delete this->config;
	delete this->world;
	delete this->random;
}

/* TODO after one run() the observables' files should 
 * be closed so that the intermediate results can be 
 * accessed (e.g. from python). is that so??? 
 * When run is called the observables' files should 
 * be opened again.*/
void Simulation::run()
{
	this->resetForces();
	this->resetActivePairs();
	world->energy = 0.;
	this->calculateInteractionForcesEnergies();
	this->calculateGeometryForcesEnergies();
	this->recordObservables(0);
	double acceptance = 1.;
	bool isStepAccepted = true;
	for (unsigned long timeIndex = 1; timeIndex < config->maxTime; timeIndex++)
	{
		/* Dynamics */
		this->saveOldState();
		this->propagateDynamics(); // propose
		this->resetForces();
		this->resetActivePairs();
		world->energy = 0.;
		this->calculateInteractionForcesEnergies(); // calculate energy and force
		this->calculateGeometryForcesEnergies();
		acceptance = this->acceptanceDynamics();
		world->acceptProbDynamics = acceptance;
		isStepAccepted = this->acceptOrReject(acceptance);
		if (config->isReversibleDynamics && ( ! isStepAccepted ) ) {
			this->restoreOldState();
			world->rejectionsDynamics += 1;
		}
		else { world->acceptionsDynamics += 1; }

		/* Reactions */
		this->saveOldState();
		acceptance = this->propagateReactions();
		this->resetForces();
		this->resetActivePairs();
		world->energy = 0.;
		this->calculateInteractionForcesEnergies();
		this->calculateGeometryForcesEnergies();
		acceptance *= this->acceptanceReactions();
		world->acceptProbReactions = acceptance;
		isStepAccepted = this->acceptOrReject(acceptance);
		if (config->isReversibleReactions && ( ! isStepAccepted ) ) {
			this->restoreOldState();
			world->rejectionsReactions += 1;
		}
		else { world->acceptionsReactions += 1; }

		/* Advance clock */
		world->cumulativeRuntime += config->timestep;
		this->recordObservables(timeIndex);
	}
}

void Simulation::saveOldState()
{
	world->oldEnergy          = world->energy;
	world->oldActiveParticles = world->activeParticles;
	world->oldActivePairs     = world->activePairs;
}

void Simulation::restoreOldState()
{
	world->energy          = world->oldEnergy;
	world->activeParticles = world->oldActiveParticles;
	world->activePairs     = world->oldActivePairs;
}

void Simulation::propagateDynamics()
{
	std::vector<double> noiseTerm = {0.,0.,0.};
	std::vector<double> forceTerm = {0.,0.,0.};
	double noisePrefactor = 1.;
	double forcePrefactor = 1.;
	double diffConst = 1.; //diffusion constant of current particle
	for (unsigned long i=0; i<world->activeParticles.size(); i++)
	{
		// look up particles' diffusion constant from its typeId
		diffConst = config->typeDict[
			world->activeParticles[i].typeId].diffusionConstant;

		noiseTerm = random->normal3D();
		noisePrefactor = sqrt(2. * diffConst * config->timestep);
		noiseTerm[0] *= noisePrefactor;
		noiseTerm[1] *= noisePrefactor;
		noiseTerm[2] *= noisePrefactor;

		forcePrefactor = config->timestep * diffConst 
			/ (config->kBoltzmann * config->temperature);
		forceTerm[0] = world->activeParticles[i].cumulativeForce[0] * forcePrefactor;
		forceTerm[1] = world->activeParticles[i].cumulativeForce[1] * forcePrefactor;
		forceTerm[2] = world->activeParticles[i].cumulativeForce[2] * forcePrefactor;

		world->activeParticles[i].move(noiseTerm);
		world->activeParticles[i].move(forceTerm);

		if (config->isPeriodic)
		{
			if (world->activeParticles[i].position[0] < (-0.5 * config->boxsize) ) {
				world->activeParticles[i].position[0] += config->boxsize;
				world->activeParticles[i].boxCoordinates[0] -= 1;
			}
			else if (world->activeParticles[i].position[0] >= (0.5 * config->boxsize) ) {
				world->activeParticles[i].position[0] -= config->boxsize;
				world->activeParticles[i].boxCoordinates[0] += 1;
			}
			if (world->activeParticles[i].position[1] < (-0.5 * config->boxsize) ) {
				world->activeParticles[i].position[1] += config->boxsize;
				world->activeParticles[i].boxCoordinates[1] -= 1;
			}
			else if (world->activeParticles[i].position[1] >= (0.5 * config->boxsize) ) {
				world->activeParticles[i].position[1] -= config->boxsize;
				world->activeParticles[i].boxCoordinates[1] += 1;
			}
			if (world->activeParticles[i].position[2] < (-0.5 * config->boxsize) ) {
				world->activeParticles[i].position[2] += config->boxsize;
				world->activeParticles[i].boxCoordinates[2] -= 1;
			}
			else if (world->activeParticles[i].position[2] >= (0.5 * config->boxsize) ) {
				world->activeParticles[i].position[2] -= config->boxsize;
				world->activeParticles[i].boxCoordinates[2] += 1;
			}
		}
	}
}

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
	for (unsigned long i=0; i<world->activePairs.size(); i++) {
		// extract indices of the pair
		unsigned long particle1Index = world->activePairs[i][0];
		unsigned long particle2Index = world->activePairs[i][1];
		// determine types of the pair
		unsigned int particle1Type
			= world->activeParticles[particle1Index].typeId;
		unsigned int particle2Type
			= world->activeParticles[particle2Index].typeId;
		std::vector<unsigned int> types = {particle1Type, particle2Type};
		for (unsigned int j=0; j<config->possibleReactions.size(); j++) {
			if ( config->possibleReactions[j]->isAffectedForward(types) ) {
				std::vector<unsigned long long> participants;
				participants.push_back(
					world->activeParticles[particle1Index].uniqueId);
				participants.push_back(
					world->activeParticles[particle2Index].uniqueId);
				ReactionEvent event(
					j, // reactionId
					true, // forwardOrBackward
					participants); // uniqueIds of reaction participants
				reactionCandidates.push_back(event);
			}
			else if ( config->possibleReactions[j]->isAffectedBackward(types) ) {
				std::vector<unsigned long long> participants;
				participants.push_back(
					world->activeParticles[particle1Index].uniqueId);
				participants.push_back(
					world->activeParticles[particle2Index].uniqueId);
				ReactionEvent event(
					j, // reactionId
					false, // forwardOrBackward
					participants); // uniqueIds of reaction participants			
				reactionCandidates.push_back(event);
			}
		}
	}
	/* Find unimolecular candidates */
	/* unimolecularParticleTypes is a list of typeIds that can undergo a
	 * unimolecular reaction like a -> b or a -> b+c, it also contains
	 * the corresponding reactionId and direction. */
	std::vector<UnimolecularCandidate> unimolecularParticleTypes;
	/* Search the possibleReactions for unimolecular reactions*/
	for (unsigned i=0; i<config->possibleReactions.size(); i++) {
		if ( config->possibleReactions[i]->forwardTypes.size() == 1 ) {
			UnimolecularCandidate uni(
				config->possibleReactions[i]->forwardTypes[0],//particleTypeId
				i, // reactionId
				true); // forwardOrBackward
			unimolecularParticleTypes.push_back(uni);
		}
		if ( config->possibleReactions[i]->backwardTypes.size() == 1 ) {
			UnimolecularCandidate uni(
				config->possibleReactions[i]->backwardTypes[0],//particleTypeId
				i, // reactionId
				false); // forwardOrBackward
			unimolecularParticleTypes.push_back(uni);
		}
	}
	/* From the above, now actually find the particular 
	 * unimolecular reaction candidates. */
	for (unsigned long i=0; i<world->activeParticles.size(); i++) {
		for (unsigned j=0; j<unimolecularParticleTypes.size(); j++) {
			/* if the type of the current particle matches a type of the
			 * prepared unimolecularParticleTypes then this unimoelcuar reaction
			 * event is added to the actual reactionCandidates. */
			if (world->activeParticles[i].typeId 
				== unimolecularParticleTypes[j].particleTypeId) {
				std::vector<unsigned long long> participants;
				participants.push_back(world->activeParticles[i].uniqueId);
				ReactionEvent event(
					unimolecularParticleTypes[j].reactionId, // reactionId
					unimolecularParticleTypes[j].forwardOrBackward,//direction
					participants); // uniqueId
				reactionCandidates.push_back(event);
			}
		}
	}
	// TODO check seeding
	std::random_shuffle(reactionCandidates.begin(), reactionCandidates.end());
	/* Perform reactions subsequently */
	double conditionalProbs = 1.;
	for (unsigned long i=0; i<reactionCandidates.size(); i++) {
		/* First find out if all participants still exist */
		bool participantsAlive = true;
		std::vector<unsigned long> particleIndices;
		std::vector<unsigned long long> participants 
			= reactionCandidates[i].participants;
		for (unsigned j=0; j<participants.size(); j++) {
			signed long index = this->findParticleIndex(participants[j]);
			if (index == -1) {participantsAlive = false;}
			else {particleIndices.push_back(index);}
		}
		/* If one of the particpants does not exist, 
		 * continue with the next reactionEvent */
		if (participantsAlive == false) {continue;}
		/* Find out the direction of reaction and perform */
		bool forwardOrBackward = reactionCandidates[i].forwardOrBackward;
		unsigned reactionId = reactionCandidates[i].reactionId;
		if (forwardOrBackward == true) {
			conditionalProbs *= config->possibleReactions[reactionId]->performForward(
				particleIndices,
				world,
				config->timestep);
		}
		if (forwardOrBackward == false) {
			conditionalProbs *= config->possibleReactions[reactionId]->performBackward(
				particleIndices,
				world,
				config->timestep);
		}
	}
	return conditionalProbs;
}

void Simulation::recordObservables(unsigned long int timeIndex)
{
	for (auto* obs : config->observables) {
		if (timeIndex % obs->recPeriod == 0) {
			obs->record(world, world->cumulativeRuntime);
		}
	}
}

void Simulation::calculateInteractionForcesEnergies()
{
	/* We want to find out if neighborlist pays off. 
	 * First determine the minimal size of subboxes,
	 * therefore check all interaction distances and
	 * all reaction radii. */
	double minimalLength = config->boxsize;
	/* check interaction distances */
	for (unsigned int i=0; i<config->possibleInteractions.size(); i++) {
		if (config->possibleInteractions[i]->cutoff < minimalLength) {
			minimalLength = config->possibleInteractions[i]->cutoff;
		}
	}
	/* check reaction radii combinations (i,j) with i <= j */
	for (unsigned int i=0; i<config->typeDict.size(); i++) {
		for (unsigned int j=i; j<config->typeDict.size(); j++) {
			double reactionRadii
				= config->typeDict[i].reactionRadius + config->typeDict[j].reactionRadius;
			if (reactionRadii < minimalLength) {
				minimalLength = reactionRadii;
			}
		}
	}
	double counter = 1.;
	unsigned int numberBoxes = 0;
	while ( (config->boxsize / counter) > minimalLength) {
		numberBoxes += 1;
		counter += 1.;
	}
	/* if n=3 we will have 9 subboxes of length L/n, which
	 * will result in having to check every box. This is
	 * as inefficient as double looping. So:
	 * ONLY construct neighborlist if we have at least 16
	 * subboxes or n>3 */
	if ( ( numberBoxes > 3 ) && config->useNeighborList ) {
		this->calculateInteractionForcesEnergiesWithLattice(numberBoxes);
	}
	else {
		this->calculateInteractionForcesEnergiesNaive();
	}
}

void Simulation::calculateInteractionForcesEnergiesNaive()
{
	for (unsigned long i=0; i<world->activeParticles.size(); i++) {
		for (unsigned long j=i+1; j<world->activeParticles.size(); j++) {
			this->calculateSingleForceEnergy(i, j);
		}
	}
}

void Simulation::calculateInteractionForcesEnergiesWithLattice(
	unsigned int numberBoxes)
{
	// construct neighborlist
	double n = (double) numberBoxes;
	double boxLength = config->boxsize / n;

	std::vector< std::vector< std::vector< std::vector<unsigned long> > > >
	neighborList(numberBoxes,
		std::vector< std::vector< std::vector<unsigned long> > >(numberBoxes,
			std::vector< std::vector<unsigned long> >(numberBoxes,
				std::vector<unsigned long>(0) ) ) );

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
	for (unsigned long j=0; j<world->activeParticles.size(); j++) {
		delX = world->activeParticles[j].position[0] + 0.5*config->boxsize;
		xIndex = (unsigned int) floor(delX / boxLength);
		delY = world->activeParticles[j].position[1] + 0.5*config->boxsize;
		yIndex = (unsigned int) floor(delY / boxLength);
		delZ = world->activeParticles[j].position[2] + 0.5*config->boxsize;
		zIndex = (unsigned int) floor(delZ / boxLength);
		// add the particles index to the list of the corresponding box
		print("push_back " << xIndex << yIndex << zIndex)
		neighborList[xIndex][yIndex][zIndex].push_back(j);
	}
	// neighborList created
	// set up vector NxN filled with bools (false initially)
	std::vector< std::vector<bool> > alreadyCalculatedPairs;
	alreadyCalculatedPairs.resize(world->activeParticles.size());
	for (unsigned int i=0; i<alreadyCalculatedPairs.size(); i++) {
		alreadyCalculatedPairs[i].resize(world->activeParticles.size());
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
		world->activeParticles[indexI].position, 
		world->activeParticles[indexJ].position, 
		config->isPeriodic, 
		config->boxsize);
	// distance of particle i,j squared
	double rSquared = r_ij[0]*r_ij[0] + r_ij[1]*r_ij[1] + r_ij[2]*r_ij[2];
	// radii of particle i and j
	double radiusI = config->typeDict[
		world->activeParticles[indexI].typeId].radius; 
	double radiusJ = config->typeDict[
		world->activeParticles[indexJ].typeId].radius;
	// squared sum of particles i,j radii
	double radiiSquared = pow(radiusI + radiusJ, 2.);
	// squared sum of reactionRadii of particle i and j
	double reactionRadiiSquared
		= pow(
		 config->typeDict[world->activeParticles[indexI].typeId].reactionRadius
		+config->typeDict[world->activeParticles[indexJ].typeId].reactionRadius
		, 2.);
	// check if particles i, j are within reactive distance
	// if so, their unique ids will be added to activePairs
	if ( rSquared <= reactionRadiiSquared ) {
		std::vector<unsigned long> activePair;
		activePair.push_back(indexI);
		activePair.push_back(indexJ);
		activePair.shrink_to_fit();
		world->activePairs.push_back(activePair);
	}
	// look for force that affects the pair (i,j)
	for (unsigned int k=0; k<config->possibleInteractions.size(); k++) {
		if (
			config->possibleInteractions[k]->isAffected(
				world->activeParticles[indexI].typeId,
				world->activeParticles[indexJ].typeId)
		) {
			// actual force call here
			config->possibleInteractions[k]->calculateForceEnergy(
				forceI,//out
				energyBuffer,
				r_ij,//in
				rSquared,
				radiiSquared);
			forceJ[0] = -1. * forceI[0];
			forceJ[1] = -1. * forceI[1];
			forceJ[2] = -1. * forceI[2];
			world->activeParticles[indexI].addForce(forceI);
			world->activeParticles[indexJ].addForce(forceJ);
			world->energy += energyBuffer;
		}
	}
}

void Simulation::calculateGeometryForcesEnergies()
{
	std::vector<double> forceI = {0.,0.,0.};
	double energyBuffer = 0.;
	for (unsigned long i=0; i<world->activeParticles.size(); i++) {
		for (unsigned long j=0; j<config->geometries.size(); j++) {
			if (config->geometries[j]->doesInteract(world->activeParticles[i].typeId)) {
				config->geometries[j]->forceEnergy(
					forceI,
					energyBuffer,
					world->activeParticles[i].position,
					config->typeDict[world->activeParticles[i].typeId].radius);
				world->activeParticles[i].addForce(forceI);
				world->energy += energyBuffer;
			}
		}
	}
}

void Simulation::resetForces()
{
	for (unsigned long i=0; i<world->activeParticles.size(); i++) {
		world->activeParticles[i].resetForce();
	}
}

void Simulation::resetActivePairs()
{
	world->activePairs.clear();
}

/* For the dynamical acceptance probability both states, old and new
 * have to have the same number of particles.  */
double Simulation::acceptanceDynamics()
{
	double acceptance = 1.;
	double firstTerm  = 0.;
	double secondTerm = 0.;
	std::vector<double> deltaX = {0.,0.,0.};
	for (unsigned long i=0; i<world->activeParticles.size(); i++) {
		getMinDistanceVector(
			deltaX,
			world->oldActiveParticles[i].position,
			world->activeParticles[i].position,
			config->isPeriodic,
			config->boxsize);
		firstTerm  += deltaX[0]
		            * ( world->oldActiveParticles[i].cumulativeForce[0]
		              + world->activeParticles[i].cumulativeForce[0] );
		firstTerm  += deltaX[1]
		            * ( world->oldActiveParticles[i].cumulativeForce[1]
		              + world->activeParticles[i].cumulativeForce[1] );
		firstTerm  += deltaX[2]
		            * ( world->oldActiveParticles[i].cumulativeForce[2]
		              + world->activeParticles[i].cumulativeForce[2] );
		secondTerm += 
		          config->typeDict[world->activeParticles[i].typeId].diffusionConstant
		            * ( world->activeParticles[i].cumulativeForce[0]
		              * world->activeParticles[i].cumulativeForce[0]
		              + world->activeParticles[i].cumulativeForce[1]
		              * world->activeParticles[i].cumulativeForce[1]
		              + world->activeParticles[i].cumulativeForce[2]
		              * world->activeParticles[i].cumulativeForce[2]
		              + world->oldActiveParticles[i].cumulativeForce[0]
		              * world->oldActiveParticles[i].cumulativeForce[0]
		              + world->oldActiveParticles[i].cumulativeForce[1]
		              * world->oldActiveParticles[i].cumulativeForce[1]
		              + world->oldActiveParticles[i].cumulativeForce[2]
		              * world->oldActiveParticles[i].cumulativeForce[2] );
	}
	firstTerm  *= 0.5;
	secondTerm *= config->timestep / (4. * config->kBoltzmann * config->temperature);
	acceptance = firstTerm + secondTerm + world->energy - world->oldEnergy;
	acceptance /= -1. * config->kBoltzmann * config->temperature;
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
	unsigned long max = world->activeParticles.size() - 1;
	unsigned long min = 0;
	unsigned long mid = 0;
	while (max >= min) {
		mid = min + (max - min) / 2;
		if (world->activeParticles[mid].uniqueId == id) {return mid;}
		else if (world->activeParticles[mid].uniqueId < id) {min = mid + 1;}
		else {max = mid - 1;}
	}
	// particle was not found
	return -1;
}