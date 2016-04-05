/* SimulationImpl.h */
#ifndef __SIMULATIONIMPL_H_INCLUDED__
#define __SIMULATIONIMPL_H_INCLUDED__
#include <math.h>
#include <cmath>
#include <vector>
#include <array>
#include <iostream>
#include <string>
#include <algorithm>
#include "World.h"
#include "Config.h"
#include "Particle.h"
#include "Random.h"
#include "Geometry.h"
#include "ParticleType.h"
#include "ParticleInteraction.h"
#include "Reaction.h"
#include "ReactionEvent.h"
#include "UnimolecularCandidateType.h"
#include "Neighborlist.h"
#include "utils.h"
#include "logging.h"
#include "Exception.h"

#include "Observable.h"
#include "Trajectory.h"
#include "TrajectoryUnique.h"
#include "RadialDistribution.h"
#include "MeanSquaredDisplacement.h"
#include "ProbabilityDensity.h"
#include "Energy.h"
#include "Acceptance.h"
#include "ParticleNumbers.h"

class SimulationImpl {
public:
	SimulationImpl(World * inWorld, Config * inConfig);
	SimulationImpl(); // default constructor if child is created
	~SimulationImpl();
	/* Start the simulation. Iterate for maxTime timesteps.*/
	void run(const unsigned long maxTime);

	/* World stores positions of particles and other variables
	 * that change during the simulation. Config stores information
	 * that does not change during the simulation. Together they
	 * contain all the system information. */
	World * world;
	Config * config;

	bool useNeighborlist;

	/* Observable methods */
	void writeAllObservablesToFile();
	std::string showObservables();
	void deleteAllObservables();
	void new_Trajectory(unsigned long recPeriod, std::string filename);
	void new_RadialDistribution(unsigned long recPeriod, std::string filename, std::vector<double> ranges, std::vector< std::array<unsigned,2> > considered);
	void new_MeanSquaredDisplacement(unsigned long recPeriod, std::string filename,	unsigned particleTypeId);
	void new_ProbabilityDensity(unsigned long recPeriod, std::string filename, unsigned pTypeId, std::vector<double> range, unsigned int coord);
	void new_Energy(unsigned long recPeriod, std::string filename);
	void new_Acceptance(unsigned long recPeriod, std::string filename, bool reactionsOrDiffusion);
	void new_ParticleNumbers(unsigned long recPeriod, std::string filename,	unsigned particleTypeId);

/* children of SimulationImpl need access to these. How to do that 
 * with private members?*/
//private: 
	/*------- core functions/variables -------*/
	Random * random; // the random number generator
	Utils * utils; // functions for calculating distance
	Neighborlist * neighborlist;
	/* If no interactions or reactions are present, no distances have to be calculated */
	bool skipPairInteractionsReactions;
	/* This assures that a Neighborlist class has been allocated. Necessary if the Simulation
	 * is destroyed before a run() has been executed. */
	bool neighborlistConfigured;
	std::vector< std::unique_ptr<Observable> > observables;
	/* unimolecularCandidateTypes is a list of typeIds that can undergo a
	 * unimolecular reaction like a -> b or a -> b+c, it also contains
	 * the corresponding reactionId and direction. */
	std::vector<UnimolecularCandidateType> unimolecularCandidateTypes;
	void configureNeighborlist();
	/* Construct the vector unimolecularCandidateTypes */
	void setupUnimolecularCandidateTypes();
	/* Configure and setup all observables, e.g. MSD gets all initial coordinates.
	 * Here configuration happens at every call on run(), whereas setup()
	 * only occurs once at the first call to run(). */
	void configureAndSetupObservables();
	/* Store the objects describing the current state:
	 * energy, particles (positions, forces) and activePairs */
	void saveOldState();
	/* Doing the opposite of the above */
	void restoreOldState();
	/* Perform Brownian Dynamics step on particles
	 * according to their accumulated forces */
	void propagateDiffusion();
	/* Perform reactions. Unimolecular according to their
	 * reaction rates and bimolecular only when the pair
	 * is part of activePairs. Already return the ratios
	 * of forward and backward probabilities. */
	double propagateReactions();
	/* Call every observables' record() function, if the timeIndex
	 * corresponds to the predefined recording interval.*/
	void recordObservables(unsigned long timeIndex);
	/* First determine how to calculate the forces, i.e. if a 
	 * neighborList approach pays off (have more than 9 boxes). */
	void calculateInteractionForcesEnergies(); 
	/* double loop (i,j) over particles and call 
	 * according Forcetype for particle pair (i,j) --> O(n^2) */
	void calculateInteractionForcesEnergiesNaive();
	/* does the same as above, but only considers interactions
	 * of neighboring boxes, that have the size of the maximum
	 * cutoff distance --> O(n) */
	void calculateInteractionForcesEnergiesWithLattice();
	/* evaluate the force and energy for a given pair of
	 * particles and store their unique ids in activePairs
	 * if they are in reactive distance */
	void calculateSingleForceEnergyCheckReactionCandidate(unsigned indexI, unsigned indexJ);
	void calculateGeometryForcesEnergies();
	void resetForces();
	void resetReactionCandidates();
	double acceptanceDiffusion();
	double acceptanceReactions();
	bool acceptOrReject(double acceptance);
	/* Return the position in particles of the 
	 * particle with the uniqueId id. Particles are added
	 * to particles only using addParticle(), which
	 * means that it is sorted w.r.t. uniqueIds. Therefore
	 * this function performs a binary search with
	 * complexity O(log n). The return value is signed so
	 * the case "no particle found" is expressed by "-1" */
	long findParticleIndex(unsigned long long id);

	/* members that are once allocated so that less allocation
	 * appears on the run */
	std::vector<double> forceI;
	std::vector<double> forceJ;
	std::vector<double> r_ij;
};
#endif//__SIMULATIONIMPL_H_INCLUDED__