/* Simulation.cpp
 * author: Christoph Froehner */

//#define __DEBUG__
#ifndef __DEBUG__
#define print(a)  
#endif
#ifdef __DEBUG__
#define print(a) std::cout << a << std::endl;
#undef __DEBUG__
#endif

#include "Simulation.h"

Simulation::Simulation()
{
	this->random = new Random("ranlxs0");
	this->world  = new World();
	this->config = new Config(this->world, this->random);
	this->utils = new Utils();
	this->forceI = {0.,0.,0.};
	this->forceJ = {0.,0.,0.};
	this->energyBuffer = 0.;
	this->r_ij = {0.,0.,0.};
	this->rSquared = 0.;
	this->radiusI = 0.;
	this->radiusJ = 0.;
	this->radiiSquared = 0.;
	this->reactionRadiiSquared = 0.;
	this->sizePossibleInteractions = 0;
}

Simulation::~Simulation()
{
	delete this->neighborlist;
	delete this->config;
	delete this->world;
	delete this->utils;
	delete this->random;
}

/* TODO after one run() the observables' files should 
 * be closed so that the intermediate results can be 
 * accessed (e.g. from python). is that so??? 
 * When run is called the observables' files should 
 * be opened again.*/
void Simulation::run()
{
	print("Enter Run")
	if (config->useNeighborList) {
		this->neighborlist = new Neighborlist(config->numberBoxes);
	}
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
		print("Enter Dynamics")
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
		print("Enter Reactions")
		this->saveOldState();
		acceptance = this->propagateReactions();
		print("After propagateReactions")
		this->resetForces();
		this->resetActivePairs();
		world->energy = 0.;
		print("After reset")
		this->calculateInteractionForcesEnergies();
		print("After calculate interactions")
		this->calculateGeometryForcesEnergies();
		print("After calculate geometries")
		acceptance *= this->acceptanceReactions();
		world->acceptProbReactions = acceptance;
		isStepAccepted = this->acceptOrReject(acceptance);
		if (config->isReversibleReactions && ( ! isStepAccepted ) ) {
			this->restoreOldState();
			world->rejectionsReactions += 1;
		}
		else { world->acceptionsReactions += 1; }

		/* Advance clock */
		print("Advance clock")
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

		if (config->isPeriodic) {
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

//#define __DEBUG__
#ifndef __DEBUG__
#define print(a)  
#endif
#ifdef __DEBUG__
#define print(a) std::cout << a << std::endl;
#undef __DEBUG__
#endif

double Simulation::propagateReactions()
{
	print("Enter propagateReactions")
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
	print("Bimolecular candidates found")
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
	print("Unimolecular types determined")
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
	print("Unimolecular candidates found")
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
			print("before find")
			signed long index = this->findParticleIndex(participants[j]);
			print("after find")
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
		else if (forwardOrBackward == false) {
			conditionalProbs *= config->possibleReactions[reactionId]->performBackward(
				particleIndices,
				world,
				config->timestep);
		}
	}
	print("Reactions performed")
	/* Check if newly created particles were put outside the box and correct */
	for (unsigned long i=0; i<world->activeParticles.size(); i++) {
		if (config->isPeriodic) {
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
	print("Checked for box-overflows")
	return conditionalProbs;
}

void Simulation::recordObservables(unsigned long timeIndex)
{
	for (unsigned i=0; i<config->observables.size(); i++) {
		if (timeIndex % config->observables[i]->recPeriod == 0) {
			config->observables[i]->record(world, world->cumulativeRuntime);
		}
	}
}

void Simulation::calculateInteractionForcesEnergies()
{
	if ( config->useNeighborList ) {
		this->calculateInteractionForcesEnergiesWithLattice(config->numberBoxes);
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
		this->neighborlist->addIndex(xIndex, yIndex, zIndex, j);
	}
	// this->neighborlist created
	// set up vector NxN filled with bools (false initially)
	/*
	std::vector< std::vector<bool> > alreadyCalculatedPairs;
	alreadyCalculatedPairs.resize(world->activeParticles.size());
	for (unsigned long i=0; i<alreadyCalculatedPairs.size(); i++) {
		alreadyCalculatedPairs[i].resize(world->activeParticles.size());
		std::fill(
			alreadyCalculatedPairs[i].begin(),
			alreadyCalculatedPairs[i].end(),
			false);
	}
	*/

	signed int otherX = 0;
	signed int otherY = 0;
	signed int otherZ = 0;
	signed x_i = 0;
	signed y_i = 0;
	signed z_i = 0;
	for (unsigned int x=0; x<numberBoxes; x++)
	for (unsigned int y=0; y<numberBoxes; y++)
	for (unsigned int z=0; z<numberBoxes; z++) {
		/* usage of half point symmetry check selfbox (x,y,z) and
		 * (x+1,y,z) and
		 * (x-1,y+1,z), (x,y+1,z), (x+1,y+1,z) and
		 * with z_i = z + 1, all 9 (x,y,z_i) pairs, with x,y in [-1,0,1] */

		//for (signed int x_i = -1; x_i < 2; x_i++)
		//for (signed int y_i = -1; y_i < 2; y_i++)
		//for (signed int z_i = -1; z_i < 2; z_i++) {
		for (unsigned k=0; k<14; k++) {
			/* these are the relative coordinates 
			 * of boxes that need to be searched
			 * in the neighborlattice method. */
			switch (k) {
				case 0: x_i= 0; y_i= 0; z_i=0;break;
				case 1: x_i= 1; y_i= 0; z_i=0;break;
				case 2: x_i=-1; y_i= 1; z_i=0;break;
				case 3: x_i= 0; y_i= 1; z_i=0;break;
				case 4: x_i= 1; y_i= 1; z_i=0;break;
				case 5: x_i=-1; y_i=-1; z_i=1;break;
				case 6: x_i= 0; y_i=-1; z_i=1;break;
				case 7: x_i= 1; y_i=-1; z_i=1;break;
				case 8: x_i=-1; y_i= 0; z_i=1;break;
				case 9: x_i= 0; y_i= 0; z_i=1;break;
				case 10:x_i= 1; y_i= 0; z_i=1;break;
				case 11:x_i=-1; y_i= 1; z_i=1;break;
				case 12:x_i= 0; y_i= 1; z_i=1;break;
				case 13:x_i= 1; y_i= 1; z_i=1;break;
			}
			if ( (x_i==0) && (y_i==0) && (z_i==0) ) {
				for (unsigned long i=0; i<this->neighborlist->getSize(x,y,z); i++)
				for (unsigned long j=i+1; j<this->neighborlist->getSize(x,y,z); j++) {
					this->calculateSingleForceEnergy(
						this->neighborlist->getIndex(x,y,z,i),
						this->neighborlist->getIndex(x,y,z,j));
					/*
					alreadyCalculatedPairs
						[this->neighborlist->getIndex(x,y,z,i)]
						[this->neighborlist->getIndex(x,y,z,j)] = true;
					alreadyCalculatedPairs
						[this->neighborlist->getIndex(x,y,z,j)]
						[this->neighborlist->getIndex(x,y,z,i)] = true;
					*/
				}
			}
			else {
				//determine the "other" subbox. detect "over/underflows"
				//TODO only detect "over/underflows" when system is periodic
				otherX = x + x_i;
				if (otherX == -1) {otherX = numberBoxes - 1;}
				if (otherX == numberBoxes) {otherX = 0;}
				otherY = y + y_i;
				if (otherY == -1) {otherY = numberBoxes - 1;}
				if (otherY == numberBoxes) {otherY = 0;}
				otherZ = z + z_i;
				if (otherZ == -1) {otherZ = numberBoxes - 1;}
				if (otherZ == numberBoxes) {otherZ = 0;}
				for (unsigned long i=0; i<this->neighborlist->getSize(x,y,z); i++)
				for (unsigned long j=0; j<this->neighborlist->getSize(otherX,otherY,otherZ); j++) {
					/*
					// if not already calculated
					if (! alreadyCalculatedPairs
						[ this->neighborlist->getIndex(x,y,z,i) ]
						[ this->neighborlist->getIndex(otherX,otherY,otherZ,j) ] ) {
						// calculate interaction
						this->calculateSingleForceEnergy(
							this->neighborlist->getIndex(x,y,z,i),
							this->neighborlist->getIndex(otherX,otherY,otherZ,j) );
						alreadyCalculatedPairs
							[ this->neighborlist->getIndex(x,y,z,i) ]
							[ this->neighborlist->getIndex(otherX,otherY,otherZ,j) ] = true;
						alreadyCalculatedPairs
							[ this->neighborlist->getIndex(otherX,otherY,otherZ,j) ]
							[ this->neighborlist->getIndex(x,y,z,i) ] = true;
					}
					*/
					// TODO verify this method
					this->calculateSingleForceEnergy(
						this->neighborlist->getIndex(x,y,z,i),
						this->neighborlist->getIndex(otherX,otherY,otherZ,j) );
				}
			}
		}
	}
	this->neighborlist->clear();
}

void Simulation::calculateSingleForceEnergy(
	unsigned int indexI,
	unsigned int indexJ)
{
	this->forceI[0]=0.; this->forceI[1]=0.; this->forceI[2]=0.;
	this->forceJ[0]=0.; this->forceJ[1]=0.; this->forceJ[2]=0.;
	// interaction energy of particle pair (i,j)
	this->energyBuffer = 0.; 
	/* First gather all parameters, distances and reaction
	 * radii of particles, check if they are within reaction
	 * distance and then check for forces that apply to the given
	 * pair of particles. */
	// connecting vector from particle i to j
	this->r_ij[0]=0.; this->r_ij[1]=0.; this->r_ij[2]=0.;
	this->utils->getMinDistanceVector(
		r_ij,
		world->activeParticles[indexI].position, 
		world->activeParticles[indexJ].position, 
		config->isPeriodic,
		config->boxsize);
	// distance of particle i,j squared
	this->rSquared = r_ij[0]*r_ij[0] + r_ij[1]*r_ij[1] + r_ij[2]*r_ij[2];
	// radii of particle i and j
	this->radiusI = config->typeDict[
		world->activeParticles[indexI].typeId].radius; 
	this->radiusJ = config->typeDict[
		world->activeParticles[indexJ].typeId].radius;
	// squared sum of particles i,j radii
	this->radiiSquared = pow(radiusI + radiusJ, 2.);
	// squared sum of reactionRadii of particle i and j
	this->reactionRadiiSquared
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
		//activePair.shrink_to_fit();
		world->activePairs.push_back(activePair);
	}
	// look for force that affects the pair (i,j)
	this->sizePossibleInteractions = config->possibleInteractions.size();
	for (unsigned int k=0; k<this->sizePossibleInteractions; k++) {
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

/* This calculates only half of the interaction between I and J
 * so that the forces/energies are not accidently calculated 
 * and added to the particles twice*/
void Simulation::calculateSingleForceEnergyOnlyForI(
	unsigned int indexI,
	unsigned int indexJ)
{
	this->forceI[0]=0.; this->forceI[1]=0.; this->forceI[2]=0.;
	// interaction energy of particle pair (i,j)
	this->energyBuffer = 0.; 
	// connecting vector from particle i to j
	this->r_ij[0]=0.; this->r_ij[1]=0.; this->r_ij[2]=0.;
	this->utils->getMinDistanceVector(
		r_ij,
		world->activeParticles[indexI].position, 
		world->activeParticles[indexJ].position, 
		config->isPeriodic, 
		config->boxsize);
	// distance of particle i,j squared
	this->rSquared = r_ij[0]*r_ij[0] + r_ij[1]*r_ij[1] + r_ij[2]*r_ij[2];
	// radii of particle i and j
	this->radiusI = config->typeDict[
		world->activeParticles[indexI].typeId].radius; 
	this->radiusJ = config->typeDict[
		world->activeParticles[indexJ].typeId].radius;
	// squared sum of particles i,j radii
	this->radiiSquared = pow(radiusI + radiusJ, 2.);
	// squared sum of reactionRadii of particle i and j
	this->reactionRadiiSquared
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
		//activePair.shrink_to_fit();
		world->activePairs.push_back(activePair);
	}
	// look for force that affects the pair (i,j)
	this->sizePossibleInteractions = config->possibleInteractions.size();
	for (unsigned int k=0; k<this->sizePossibleInteractions; k++) {
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
			world->activeParticles[indexI].addForce(forceI);
			world->energy += 0.5*energyBuffer;
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
		this->utils->getMinDistanceVector(
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
		              - world->oldActiveParticles[i].cumulativeForce[0]
		              * world->oldActiveParticles[i].cumulativeForce[0]
		              - world->oldActiveParticles[i].cumulativeForce[1]
		              * world->oldActiveParticles[i].cumulativeForce[1]
		              - world->oldActiveParticles[i].cumulativeForce[2]
		              * world->oldActiveParticles[i].cumulativeForce[2] );
	}
	firstTerm  *= 0.5;
	secondTerm *= config->timestep / (4. * config->kBoltzmann * config->temperature);
	acceptance = firstTerm + secondTerm + world->energy - world->oldEnergy;
	acceptance /= -1. * config->kBoltzmann * config->temperature;
	acceptance = exp( acceptance );
	return acceptance;
}

double Simulation::acceptanceReactions()
{
	double unconditional = 1.;
	unconditional = world->energy - world->oldEnergy;
	unconditional /= -1. * config->kBoltzmann * config->temperature;
	unconditional = exp( unconditional );
	return unconditional;
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

//#define __DEBUG__
#ifndef __DEBUG__
#define print(a)  
#endif
#ifdef __DEBUG__
#define print(a) std::cout << a << std::endl;
#undef __DEBUG__
#endif

long int Simulation::findParticleIndex(unsigned long long id)
{
	print("Enter findParticleIndex")
	signed long max = world->activeParticles.size() - 1;
	signed long min = 0;
	signed long mid = 0;
	while (max >= min) {
		mid = min + (max - min) / 2;
		print("mid set. min, mid, max: " << min << ", " << mid << ", " << max)
		if (world->activeParticles[mid].uniqueId == id) {
			return mid;
		}
		else if (world->activeParticles[mid].uniqueId < id) {min = mid + 1;}
		else {max = mid - 1;}
	}
	// particle was not found
	return -1;
}