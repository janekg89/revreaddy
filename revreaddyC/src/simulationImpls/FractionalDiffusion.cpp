#include "FractionalDiffusion.h"
#include "generateIncrements.h"

FractionalDiffusion::FractionalDiffusion(World * inWorld, Config * inConfig) {
	LOG_TRACE("Enter FractionalDiffusion constructor")
	this->world  = inWorld;
	this->config = inConfig;
	this->random = new Random("ranlxs0");
	this->utils = new Utils();
	this->forceI = {0.,0.,0.};
	this->forceJ = {0.,0.,0.};
	this->r_ij = {0.,0.,0.};
	this->useNeighborlist = true;
	this->skipPairInteractionsReactions = false;
	LOG_TRACE("Leave FractionalDiffusion constructor")
}

void FractionalDiffusion::run(const unsigned long maxTime) {
	config->maxTime = maxTime;
	config->configureReactions();
	this->skipPairInteractionsReactions = false;
	if (config->interactions.empty() && config->reactions.empty()) {
		this->skipPairInteractionsReactions = true;
	}
	if (this->useNeighborlist && ( !this->skipPairInteractionsReactions ) ) { this->configureNeighborlist(); }
	else { this->useNeighborlistThisRun = false; }
	this->setupUnimolecularCandidateTypes();
	this->configureAndSetupObservables();

	/* Set up increments */
	for (auto&& p : world->particles) {
		double diffConst = config->particleTypes[p.typeId].diffusionConstant;
		double timestep = config->timestep;
		double alpha = 0.5;
		world->increments[p.uniqueId] = janek::generateIncrements(maxTime, diffConst,timestep, alpha , random);//generateIncrements(args, this->random);
	    world->incrementsIndex[p.uniqueId] = 0;
	}

	this->resetForces();
	this->resetReactionCandidates();
	world->energy = 0.;
	this->calculateInteractionForcesEnergies();
	this->calculateGeometryForcesEnergies();
	this->recordObservables(0);
	for (unsigned long timeIndex = 0; timeIndex < maxTime; ++timeIndex) {
		/* Diffusion */
		this->propagateDiffusion(); // propose
		this->resetForces();
		this->resetReactionCandidates();
		world->energy = 0.;
		if ( !this->skipPairInteractionsReactions ) {
			this->calculateInteractionForcesEnergies(); // calculate energy and force
		}
		this->calculateGeometryForcesEnergies();
		/* Reactions */
		this->propagateReactions();
		this->resetForces();
		this->resetReactionCandidates();
		world->energy = 0.;
		if ( !this->skipPairInteractionsReactions ) {
			this->calculateInteractionForcesEnergies();
		}
		this->calculateGeometryForcesEnergies();

		/* Advance clock */
		world->cumulativeRuntime += config->timestep;
		/* +1 because recordObservables(0) was called already before the run.
		 * This means, if you call run(2): the observables will have 3 entries
		 * one for the initial state and two for the states after executing two
		 * timesteps */
		this->recordObservables(timeIndex + 1);
	}
	// clean up after run
    unimolecularCandidateTypes.clear();
    if (this->useNeighborlistThisRun) {
        delete this->neighborlist;
    }
}

void FractionalDiffusion::propagateDiffusion() {
	std::vector<double> noiseTerm = {0.,0.,0.};
	std::vector<double> forceTerm = {0.,0.,0.};
	//double noisePrefactor;
	double forcePrefactor;
	double diffConst; //diffusion constant of current particle
	for (unsigned long i=0; i<world->particles.size(); i++)
	{
		// look up particles' diffusion constant from its typeId
		diffConst = config->particleTypes[world->particles[i].typeId].diffusionConstant;

        /* Change this */
		boost::multi_array<double,2> incrementsParticle = world->increments.find(i)->second;

        long int timeIndex = world->incrementsIndex.find(i)->second;
        boost::array<boost::multi_array<double, 2> ::index,2> idx1 = {{0,timeIndex}};
        boost::array<boost::multi_array<double, 2> ::index,2> idx2 = {{1,timeIndex}};
        boost::array<boost::multi_array<double, 2> ::index,2> idx3 = {{2,timeIndex}};



        noiseTerm[0]= incrementsParticle(idx1);
        noiseTerm[1]= incrementsParticle(idx2);
        noiseTerm[2]= incrementsParticle(idx3);
        /*
		noiseTerm[0] = incrementsParticle[0][timeIndex];
		noiseTerm[1] = incrementsParticle[1][timeIndex];
		noiseTerm[2] = incrementsParticle[2][timeIndex ];
        */
        world->incrementsIndex.find(i)->second = (long unsigned int)++timeIndex;



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
