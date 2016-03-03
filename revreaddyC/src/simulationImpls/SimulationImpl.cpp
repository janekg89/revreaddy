/* SimulationImpl.cpp
 * author: Christoph Froehner */

#include "SimulationImpl.h"

template<typename T, typename ...Args>
std::unique_ptr<T> make_unique( Args&& ...args ) {
	return std::unique_ptr<T>( new T( std::forward<Args>(args)... ) );
}

SimulationImpl::SimulationImpl(World * inWorld, Config * inConfig) {
	LOG_TRACE("Enter SimulationImpl Constructor.")
	this->world  = inWorld;
	this->config = inConfig;
	this->random = new Random("ranlxs0");
	this->utils = new Utils();
	this->forceI = {0.,0.,0.};
	this->forceJ = {0.,0.,0.};
	this->r_ij = {0.,0.,0.};
	this->useNeighborlist = true;
	this->neighborlistConfigured = false;
	this->skipPairInteractionsReactions = false;
	LOG_TRACE("Leave SimulationImpl Constructor.")
}

// default constructor if child is created
SimulationImpl::SimulationImpl() {
	LOG_TRACE("Enter and leave SimulationImpl default constructor");
}

SimulationImpl::~SimulationImpl() {
	LOG_TRACE("Enter SimulationImpl Destructor.")
	this->deleteAllObservables();
	if (this->neighborlistConfigured) {
		delete this->neighborlist;
	}
	delete this->utils;
	delete this->random;
	LOG_TRACE("Leave SimulationImpl Destructor.")
}

void SimulationImpl::writeAllObservablesToFile() {
	for (unsigned i=0; i<this->observables.size(); ++i) {
		this->observables[i]->writeBufferToFile();
	}
}

void SimulationImpl::writeLastObservableToFile() {
	if (this->observables.size() > 0) {
		this->observables.back()->writeBufferToFile();
	}
	else {
		throw Exception("There are no observables to write");
	}
}

std::string SimulationImpl::showObservables() {
	std::string content = "Observables: ";
	if (this->observables.size() > 0) {
		for (unsigned i=0; i<this->observables.size(); ++i) {
			content += std::string (typeid(* this->observables[i]).name()) + " ";
		}
	}
	else {content += "empty";}
	content += "\n";
	return content;
}

void SimulationImpl::deleteAllObservables() {
	/* Erase all the unique pointers, the observables are thus deleted as well */
	this->observables.clear();
}

void SimulationImpl::deleteLastObservable() {
	/* delete the smart pointer, upon which the observable will be destroyed as well */
	if (this->observables.size() > 0) {
		this->observables.pop_back();
	}
}

void SimulationImpl::new_Trajectory(unsigned long recPeriod, std::string filename) {
	std::unique_ptr<Trajectory> obs = make_unique<Trajectory>(
		recPeriod,
		0,
		filename);
	this->observables.push_back( std::move(obs) );
}

void SimulationImpl::new_RadialDistribution(unsigned long recPeriod, std::string filename, std::vector<double> ranges, std::vector< std::vector<unsigned> > considered) {
	std::unique_ptr<RadialDistribution> rad = make_unique<RadialDistribution>(
		recPeriod,
		0,
		ranges,
		considered,
		filename);
	this->observables.push_back( std::move(rad) );
}

void SimulationImpl::new_MeanSquaredDisplacement(unsigned long recPeriod, std::string filename,	unsigned particleTypeId) {
	std::unique_ptr<MeanSquaredDisplacement> msd = make_unique<MeanSquaredDisplacement>(
		recPeriod,
		0,
		particleTypeId,
		filename);
	this->observables.push_back( std::move(msd) );
}

void SimulationImpl::new_ProbabilityDensity(unsigned long recPeriod, std::string filename, unsigned particleTypeId, std::vector<double> range, unsigned coord) {
	std::unique_ptr<ProbabilityDensity> prob = make_unique<ProbabilityDensity>(
		recPeriod,
		0,
		filename,
		particleTypeId,
		range,
		coord);
	this->observables.push_back( std::move(prob) );
}

void SimulationImpl::new_Energy(unsigned long recPeriod, std::string filename) {
	std::unique_ptr<Energy> ener = make_unique<Energy>(
		recPeriod,
		0,
		filename);
	this->observables.push_back( std::move(ener) );
}

void SimulationImpl::new_Acceptance(unsigned long recPeriod, std::string filename, bool reactionsOrDynamics) {
	std::unique_ptr<Acceptance> acc = make_unique<Acceptance>(
		recPeriod,
		0,
		filename,
		reactionsOrDynamics);
	this->observables.push_back( std::move(acc) );
}

void SimulationImpl::new_ParticleNumbers(unsigned long recPeriod, std::string filename,	unsigned particleTypeId) {
	std::unique_ptr<ParticleNumbers> par = make_unique<ParticleNumbers>(
		recPeriod,
		0,
		filename,
		particleTypeId);
	this->observables.push_back( std::move(par) );
}

void SimulationImpl::run(const unsigned long maxTime) {
	config->configureReactions();
	if (config->interactions.empty() && config->reactions.empty()) {
		this->skipPairInteractionsReactions = true;
	}
	if (this->useNeighborlist && ( !this->skipPairInteractionsReactions ) ) { this->configureNeighborlist(); }
	else { this->useNeighborlist = false; }
	this->setupUnimolecularCandidateTypes();
	this->configureObservables();
	this->resetForces();
	this->resetReactionCandidates();
	world->energy = 0.;
	this->calculateInteractionForcesEnergies();
	this->calculateGeometryForcesEnergies();
	this->recordObservables(0);
	double acceptance = 1.;
	bool isStepAccepted = true;
	for (unsigned long timeIndex = 1; timeIndex < maxTime; ++timeIndex) {
		/* Dynamics */
		this->saveOldState();
		this->propagateDynamics(); // propose
		this->resetForces();
		this->resetReactionCandidates();
		world->energy = 0.;
		if ( !this->skipPairInteractionsReactions ) {
			this->calculateInteractionForcesEnergies(); // calculate energy and force
		}
		this->calculateGeometryForcesEnergies();
		acceptance = this->acceptanceDynamics();
		world->acceptProbDynamics = acceptance;
		isStepAccepted = this->acceptOrReject(acceptance);
		// TODO revDynamics
		if (true && ( ! isStepAccepted ) ) {
			this->restoreOldState();
			world->rejectionsDynamics += 1;
		}
		else { world->acceptionsDynamics += 1; }

		/* Reactions */
		this->saveOldState();
		acceptance = this->propagateReactions();
		this->resetForces();
		this->resetReactionCandidates();
		world->energy = 0.;
		if ( !this->skipPairInteractionsReactions ) {
			this->calculateInteractionForcesEnergies();
		}
		this->calculateGeometryForcesEnergies();
		acceptance *= this->acceptanceReactions();
		world->acceptProbReactions = acceptance;
		isStepAccepted = this->acceptOrReject(acceptance);
		// TODO revReactions
		if (true && ( ! isStepAccepted ) ) {
			this->restoreOldState();
			world->rejectionsReactions += 1;
		}
		else { world->acceptionsReactions += 1; }

		/* Advance clock */
		world->cumulativeRuntime += config->timestep;
		this->recordObservables(timeIndex);
	}
}

void SimulationImpl::configureNeighborlist() {
	double minimalLength = config->boxsize;
	// check interaction distances
	for (unsigned i=0; i<config->interactions.size(); ++i) {
		if (config->interactions[i]->cutoff < minimalLength) {
			minimalLength = config->interactions[i]->cutoff;
		}
	}
	// check reaction distances
	for (unsigned i=0; i<config->reactions.size(); ++i) {
		if (config->reactions[i]->reactionDistance < minimalLength) {
			minimalLength = config->reactions[i]->reactionDistance;
		}
	}
	double counter = 1.;
	unsigned numberBoxes = 0;
	while ((config->boxsize / counter) > minimalLength) {
		numberBoxes += 1;
		counter += 1.;
	}
	if (numberBoxes > 3) {
		this->useNeighborlist = true; // actually obsolete
		this->neighborlist = new Neighborlist( numberBoxes );
		this->neighborlistConfigured = true;
	}
	else {
		this->useNeighborlist = false;
	}	
}

void SimulationImpl::setupUnimolecularCandidateTypes() {
	/* Search the reactions for unimolecular reactions */
	for (unsigned i=0; i<config->reactions.size(); i++) {
		if ( config->reactions[i]->forwardTypes.size() == 1 ) {
			UnimolecularCandidateType uni(
				config->reactions[i]->forwardTypes[0],//particleTypeId
				i, // reactionId
				true); // forwardOrBackward
			unimolecularCandidateTypes.push_back(uni);
		}
		if ( config->reactions[i]->backwardTypes.size() == 1 ) {
			UnimolecularCandidateType uni(
				config->reactions[i]->backwardTypes[0],//particleTypeId
				i, // reactionId
				false); // forwardOrBackward
			unimolecularCandidateTypes.push_back(uni);
		}
	}	
}

void SimulationImpl::saveOldState() {
	world->oldEnergy = world->energy;
	world->oldParticles = world->particles;
	world->oldReactionCandidates = world->reactionCandidates;
}

void SimulationImpl::restoreOldState() {
	world->energy = world->oldEnergy;
	world->particles = world->oldParticles;
	world->reactionCandidates = world->oldReactionCandidates;
}

void SimulationImpl::propagateDynamics() {
	std::vector<double> noiseTerm = {0.,0.,0.};
	std::vector<double> forceTerm = {0.,0.,0.};
	double noisePrefactor = 1.;
	double forcePrefactor = 1.;
	double diffConst = 1.; //diffusion constant of current particle
	for (unsigned long i=0; i<world->particles.size(); i++)
	{
		// look up particles' diffusion constant from its typeId
		diffConst = config->particleTypes[world->particles[i].typeId].diffusionConstant;

		noiseTerm = random->normal3D();
		noisePrefactor = sqrt(2. * diffConst * config->timestep);
		noiseTerm[0] *= noisePrefactor;
		noiseTerm[1] *= noisePrefactor;
		noiseTerm[2] *= noisePrefactor;

		forcePrefactor = config->timestep * diffConst / config->kT;
		forceTerm[0] = world->particles[i].cumulativeForce[0] * forcePrefactor;
		forceTerm[1] = world->particles[i].cumulativeForce[1] * forcePrefactor;
		forceTerm[2] = world->particles[i].cumulativeForce[2] * forcePrefactor;

		world->particles[i].move(noiseTerm);
		world->particles[i].move(forceTerm);

		if (config->isPeriodic) {
			if (world->particles[i].position[0] < (-0.5 * config->boxsize) ) {
				world->particles[i].position[0] += config->boxsize;
				world->particles[i].boxCoordinates[0] -= 1;
			}
			else if (world->particles[i].position[0] >= (0.5 * config->boxsize) ) {
				world->particles[i].position[0] -= config->boxsize;
				world->particles[i].boxCoordinates[0] += 1;
			}
			if (world->particles[i].position[1] < (-0.5 * config->boxsize) ) {
				world->particles[i].position[1] += config->boxsize;
				world->particles[i].boxCoordinates[1] -= 1;
			}
			else if (world->particles[i].position[1] >= (0.5 * config->boxsize) ) {
				world->particles[i].position[1] -= config->boxsize;
				world->particles[i].boxCoordinates[1] += 1;
			}
			if (world->particles[i].position[2] < (-0.5 * config->boxsize) ) {
				world->particles[i].position[2] += config->boxsize;
				world->particles[i].boxCoordinates[2] -= 1;
			}
			else if (world->particles[i].position[2] >= (0.5 * config->boxsize) ) {
				world->particles[i].position[2] -= config->boxsize;
				world->particles[i].boxCoordinates[2] += 1;
			}
		}
	}
}

double SimulationImpl::propagateReactions() {
	/* Find unimolecular candidates because bimolecular were already
	 * determined when calculating forces. */
	for (unsigned long i=0; i<world->particles.size(); i++) {
		for (unsigned j=0; j<unimolecularCandidateTypes.size(); j++) {
			/* if the type of the current particle matches a type of the
			 * prepared unimolecularCandidateTypes then this unimoelcuar reaction
			 * event is added to the actual reactionCandidates. */
			if (world->particles[i].typeId == unimolecularCandidateTypes[j].particleTypeId) {
				std::vector<unsigned long long> participants;
				participants.push_back(world->particles[i].uniqueId);
				ReactionEvent event(
					unimolecularCandidateTypes[j].reactionId, // reactionId
					unimolecularCandidateTypes[j].forwardOrBackward,//direction
					participants); // uniqueId
				world->reactionCandidates.push_back(event);
			}
		}
	}
	/* Shuffle the list reactionCandidates and perform reactions.
	 * If one particle
	 * cannot be found via its uniqueId, because it was destroyed in
	 * a previous reaction, skip this event and check the next 
	 * candidate. Note that a ReactionEvent will be performed with
	 * its predefined probability. */
	// TODO check seeding
	std::random_shuffle(world->reactionCandidates.begin(), world->reactionCandidates.end());
	/* Perform reactions subsequently */
	double conditionalProbs = 1.;
	for (unsigned long i=0; i < world->reactionCandidates.size(); i++) {
		/* First find out if all participants still exist */
		bool participantsAlive = true;
		std::vector<unsigned long> particleIndices;
		std::vector<unsigned long long> participants = world->reactionCandidates[i].participants;
		for (unsigned j=0; j<participants.size(); j++) {
			signed long index = this->findParticleIndex(participants[j]);
			if (index == -1) {participantsAlive = false;}
			else {particleIndices.push_back(index);}
		}
		/* If one of the particpants does not exist, 
		 * continue with the next reactionEvent */
		if (participantsAlive == false) {continue;}
		/* Find out the direction of reaction and perform */
		bool forwardOrBackward = world->reactionCandidates[i].forwardOrBackward;
		unsigned reactionId = world->reactionCandidates[i].reactionId;
		if (forwardOrBackward == true) {
			conditionalProbs *= config->reactions[reactionId]->performForward(
				particleIndices,
				config->timestep,
				world,
				this->random);
		}
		else if (forwardOrBackward == false) {
			conditionalProbs *= config->reactions[reactionId]->performBackward(
				particleIndices,
				config->timestep,
				world,
				this->random);
		}
	}
	/* Check if newly created particles were put outside the box and correct */
	for (unsigned long i=0; i<world->particles.size(); i++) {
		if (config->isPeriodic) {
			if (world->particles[i].position[0] < (-0.5 * config->boxsize) ) {
				world->particles[i].position[0] += config->boxsize;
				world->particles[i].boxCoordinates[0] -= 1;
			}
			else if (world->particles[i].position[0] >= (0.5 * config->boxsize) ) {
				world->particles[i].position[0] -= config->boxsize;
				world->particles[i].boxCoordinates[0] += 1;
			}
			if (world->particles[i].position[1] < (-0.5 * config->boxsize) ) {
				world->particles[i].position[1] += config->boxsize;
				world->particles[i].boxCoordinates[1] -= 1;
			}
			else if (world->particles[i].position[1] >= (0.5 * config->boxsize) ) {
				world->particles[i].position[1] -= config->boxsize;
				world->particles[i].boxCoordinates[1] += 1;
			}
			if (world->particles[i].position[2] < (-0.5 * config->boxsize) ) {
				world->particles[i].position[2] += config->boxsize;
				world->particles[i].boxCoordinates[2] -= 1;
			}
			else if (world->particles[i].position[2] >= (0.5 * config->boxsize) ) {
				world->particles[i].position[2] -= config->boxsize;
				world->particles[i].boxCoordinates[2] += 1;
			}
		}
	}
	return conditionalProbs;
}

void SimulationImpl::recordObservables(unsigned long timeIndex) {
	for (unsigned i=0; i<this->observables.size(); i++) {
		if (timeIndex % this->observables[i]->recPeriod == 0) {
			this->observables[i]->record(world, world->cumulativeRuntime);
		}
	}
}

void SimulationImpl::calculateInteractionForcesEnergies() {
	if ( this->useNeighborlist ) {
		this->calculateInteractionForcesEnergiesWithLattice();
	}
	else {
		this->calculateInteractionForcesEnergiesNaive();
	}
}

void SimulationImpl::calculateInteractionForcesEnergiesNaive() {
	for (unsigned long i=0; i<world->particles.size(); i++) {
		for (unsigned long j=i+1; j<world->particles.size(); j++) {
			this->calculateSingleForceEnergyCheckReactionCandidate(i, j);
		}
	}
}

void SimulationImpl::calculateInteractionForcesEnergiesWithLattice() {
	unsigned numberBoxes = this->neighborlist->numberBoxes;
	double n = (double) numberBoxes;
	double boxLength = config->boxsize / n;

	double delX = 0.;
	double delY = 0.;
	double delZ = 0.;
	unsigned int xIndex = 0;
	unsigned int yIndex = 0;
	unsigned int zIndex = 0;
	// find the right box triplet [x][y][z] for each particle
	for (unsigned long j=0; j<world->particles.size(); j++) {
		delX = world->particles[j].position[0] + 0.5*config->boxsize;
		xIndex = (unsigned int) floor(delX / boxLength);
		delY = world->particles[j].position[1] + 0.5*config->boxsize;
		yIndex = (unsigned int) floor(delY / boxLength);
		delZ = world->particles[j].position[2] + 0.5*config->boxsize;
		zIndex = (unsigned int) floor(delZ / boxLength);
		// add the particles index to the list of the corresponding box
		this->neighborlist->addIndex(xIndex, yIndex, zIndex, j);
	}
	// this->neighborlist created

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
					this->calculateSingleForceEnergyCheckReactionCandidate(
						this->neighborlist->getIndex(x,y,z,i),
						this->neighborlist->getIndex(x,y,z,j));
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
					this->calculateSingleForceEnergyCheckReactionCandidate(
						this->neighborlist->getIndex(x,y,z,i),
						this->neighborlist->getIndex(otherX,otherY,otherZ,j) );
				}
			}
		}
	}
	this->neighborlist->clear();
}

void SimulationImpl::calculateSingleForceEnergyCheckReactionCandidate(unsigned int indexI, unsigned int indexJ) {
	this->forceI[0]=0.; this->forceI[1]=0.; this->forceI[2]=0.;
	this->forceJ[0]=0.; this->forceJ[1]=0.; this->forceJ[2]=0.;
	// interaction energy of particle pair (i,j)
	double energyBuffer = 0.; 
	/* First gather all parameters, distances and reaction
	 * radii of particles, check if they are within reaction
	 * distance and then check for forces that apply to the given
	 * pair of particles. */
	// connecting vector from particle i to j
	this->r_ij[0]=0.; this->r_ij[1]=0.; this->r_ij[2]=0.;
	this->utils->getMinDistanceVector(
		r_ij,
		world->particles[indexI].position, 
		world->particles[indexJ].position, 
		config->isPeriodic,
		config->boxsize);
	// distance of particle i,j squared
	double rSquared = r_ij[0]*r_ij[0] + r_ij[1]*r_ij[1] + r_ij[2]*r_ij[2];
	// types of particles
	unsigned typeI = world->particles[indexI].typeId;
	unsigned typeJ = world->particles[indexJ].typeId;
	// radii of particle i and j
	double radiusI = config->particleTypes[typeI].radius; 
	double radiusJ = config->particleTypes[typeJ].radius;
	// squared sum of particles i,j radii
	double radiiSquared = pow(radiusI + radiusJ, 2.);
	// check if particles i, j are within reactive distance
	// if so, a ReactionEvent will be added to reactionCandidates
	// Therefore check for reaction that affects the pair (i,j)
	std::vector<unsigned int> types = {typeI, typeJ};
	for (unsigned k=0; k<config->reactions.size(); ++k) {
		if ( config->reactions[k]->isAffectedForward(types) ) {
			std::vector<unsigned long long> participants;
			participants.push_back(world->particles[indexI].uniqueId);
			participants.push_back(world->particles[indexJ].uniqueId);
			ReactionEvent event(
				k, // reactionId
				true, // forwardOrBackward
				participants); // uniqueIds of reaction participants
			world->reactionCandidates.push_back(event);
		}
		else if ( config->reactions[k]->isAffectedBackward(types) ) {
			std::vector<unsigned long long> participants;
			participants.push_back(world->particles[indexI].uniqueId);
			participants.push_back(world->particles[indexJ].uniqueId);
			ReactionEvent event(
				k, // reactionId
				false, // forwardOrBackward
				participants); // uniqueIds of reaction participants			
			world->reactionCandidates.push_back(event);
		}
	}
	// look for force that affects the pair (i,j)
	for (unsigned int k=0; k<config->interactions.size(); ++k) {
		if (config->interactions[k]->isAffected(
			world->particles[indexI].typeId,
			world->particles[indexJ].typeId) 
		) {
			// actual force call here
			config->interactions[k]->calculateForceEnergy(
				forceI,//out
				energyBuffer,
				r_ij,//in
				rSquared,
				radiiSquared);
			forceJ[0] = -1. * forceI[0];
			forceJ[1] = -1. * forceI[1];
			forceJ[2] = -1. * forceI[2];
			world->particles[indexI].addForce(forceI);
			world->particles[indexJ].addForce(forceJ);
			world->energy += energyBuffer;
		}
	}
}

void SimulationImpl::calculateGeometryForcesEnergies() {
	std::vector<double> forceI = {0.,0.,0.};
	double energyBuffer = 0.;
	for (unsigned long i=0; i<world->particles.size(); i++) {
		for (unsigned long j=0; j<config->geometries.size(); j++) {
			if (config->geometries[j]->doesInteract(world->particles[i].typeId)) {
				config->geometries[j]->forceEnergy(
					forceI,
					energyBuffer,
					world->particles[i].position,
					config->particleTypes[world->particles[i].typeId].radius);
				world->particles[i].addForce(forceI);
				world->energy += energyBuffer;
			}
		}
	}
}

void SimulationImpl::resetForces() {
	for (unsigned long i=0; i<world->particles.size(); i++) {
		world->particles[i].resetForce();
	}
}

void SimulationImpl::resetReactionCandidates() {
	world->reactionCandidates.clear();
}

/* For the dynamical acceptance probability both states, old and new
 * have to have the same number of particles.  */
double SimulationImpl::acceptanceDynamics() {
	double acceptance = 1.;
	double firstTerm  = 0.;
	double secondTerm = 0.;
	std::vector<double> deltaX = {0.,0.,0.};
	for (unsigned long i=0; i<world->particles.size(); i++) {
		this->utils->getMinDistanceVector(
			deltaX,
			world->oldParticles[i].position,
			world->particles[i].position,
			config->isPeriodic,
			config->boxsize);
		firstTerm  += deltaX[0]
		            * ( world->oldParticles[i].cumulativeForce[0]
		              + world->particles[i].cumulativeForce[0] );
		firstTerm  += deltaX[1]
		            * ( world->oldParticles[i].cumulativeForce[1]
		              + world->particles[i].cumulativeForce[1] );
		firstTerm  += deltaX[2]
		            * ( world->oldParticles[i].cumulativeForce[2]
		              + world->particles[i].cumulativeForce[2] );
		secondTerm += 
		          config->particleTypes[world->particles[i].typeId].diffusionConstant
		            * ( world->particles[i].cumulativeForce[0]
		              * world->particles[i].cumulativeForce[0]
		              + world->particles[i].cumulativeForce[1]
		              * world->particles[i].cumulativeForce[1]
		              + world->particles[i].cumulativeForce[2]
		              * world->particles[i].cumulativeForce[2]
		              - world->oldParticles[i].cumulativeForce[0]
		              * world->oldParticles[i].cumulativeForce[0]
		              - world->oldParticles[i].cumulativeForce[1]
		              * world->oldParticles[i].cumulativeForce[1]
		              - world->oldParticles[i].cumulativeForce[2]
		              * world->oldParticles[i].cumulativeForce[2] );
	}
	firstTerm  *= 0.5;
	secondTerm *= config->timestep / (4. * config->kT);
	acceptance = firstTerm + secondTerm + world->energy - world->oldEnergy;
	acceptance /= -1. * config->kT;
	acceptance = exp( acceptance );
	return acceptance;
}

double SimulationImpl::acceptanceReactions() {
	double unconditional = 1.;
	unconditional = world->energy - world->oldEnergy;
	unconditional /= -1. * config->kT;
	unconditional = exp( unconditional );
	return unconditional;
}

bool SimulationImpl::acceptOrReject(double acceptance) {
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

long int SimulationImpl::findParticleIndex(unsigned long long id) {
	// max min mid have to be signed, otherwise max might underflow undetected
	signed long max = world->particles.size() - 1;
	signed long min = 0;
	signed long mid = 0;
	while (max >= min) {
		mid = min + (max - min) / 2;
		if (world->particles[mid].uniqueId == id) {
			return mid;
		}
		else if (world->particles[mid].uniqueId < id) {min = mid + 1;}
		else {max = mid - 1;}
	}
	// particle was not found
	return -1;
}

void SimulationImpl::configureObservables() {
	for (unsigned i=0; i<this->observables.size(); ++i) {
		this->observables[i]->configure(world, config);
	}
}