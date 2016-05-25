#include "RevDiffusion.h"

RevDiffusion::RevDiffusion(World * inWorld, Config * inConfig) {
	LOG_TRACE("Enter RevDiffusion constructor")
	this->world  = inWorld;
	this->config = inConfig;
	this->random = new Random("ranlxs0");
	this->utils = new Utils();
	this->forceI = {0.,0.,0.};
	this->forceJ = {0.,0.,0.};
	this->r_ij = {0.,0.,0.};
	this->useNeighborlist = true;
	this->skipPairInteractionsReactions = false;
	LOG_TRACE("Leave RevDiffusion constructor")
}

void RevDiffusion::run(const unsigned long maxTime) {
	LOG_INFO("Start run() of RevDiffusion implementation.");
	config->configureReactions();
	this->skipPairInteractionsReactions = false;
	if (config->interactions.empty() && config->reactions.empty()) {
		this->skipPairInteractionsReactions = true;
	}
	if (this->useNeighborlist && ( !this->skipPairInteractionsReactions ) ) { this->configureNeighborlist(); }
	else { this->useNeighborlistThisRun = false; }
	this->setupUnimolecularCandidateTypes();
	this->configureAndSetupObservables();
	this->resetForces();
	this->resetReactionCandidates();
	world->energy = 0.;
	this->calculateInteractionForcesEnergies();
	this->calculateGeometryForcesEnergies();
	this->recordObservables(0);
	double acceptance = 1.;
	bool isStepAccepted = true;
	for (unsigned long timeIndex = 1; timeIndex < maxTime; ++timeIndex) {
		/* Diffusion */
		this->saveOldState();
		this->propagateDiffusion(); // propose
		this->resetForces();
		this->resetReactionCandidates();
		world->energy = 0.;
		if ( !this->skipPairInteractionsReactions ) {
			this->calculateInteractionForcesEnergies(); // calculate energy and force
		}
		this->calculateGeometryForcesEnergies();
		acceptance = this->acceptanceDiffusion();
		world->acceptProbDiffusion = acceptance;
		isStepAccepted = this->acceptOrReject(acceptance);
		if ( ! isStepAccepted ) {
			this->restoreOldState();
			world->rejectionsDiffusion += 1;
		}
		else { world->acceptionsDiffusion += 1; }

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
		// no acceptance step for reactions here

		/* Advance clock */
		world->cumulativeRuntime += config->timestep;
		this->recordObservables(timeIndex + 1);
	}
	// clean up after run
	unimolecularCandidateTypes.clear();
    if (this->useNeighborlistThisRun) {
        delete this->neighborlist;
    }
}