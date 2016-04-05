#include "RevReactions.h"

RevReactions::RevReactions(World * inWorld, Config * inConfig) {
	LOG_TRACE("Enter RevReactions constructor")
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
	LOG_TRACE("Leave RevReactions constructor")
}

void RevReactions::run(const unsigned long maxTime) {
	LOG_INFO("Start run() of RevReactions implementation.")
	config->configureReactions();
	if (config->interactions.empty() && config->reactions.empty()) {
		this->skipPairInteractionsReactions = true;
	}
	if (this->useNeighborlist && ( !this->skipPairInteractionsReactions ) ) { this->configureNeighborlist(); }
	else { this->useNeighborlist = false; }
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
		// no acceptance step for Diffusion here

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
		if ( ! isStepAccepted ) {
			this->restoreOldState();
			world->rejectionsReactions += 1;
		}
		else { world->acceptionsReactions += 1; }

		/* Advance clock */
		world->cumulativeRuntime += config->timestep;
		this->recordObservables(timeIndex);
	}
}