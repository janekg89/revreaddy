/* Simulation.h
 * author: Christoph Froehner
 *
 * This is the main class, which uses all other classes.
 * Simulation performs the time loop which propagates
 * the particles and also holds information that concern
 * a particular realization, e.g. the list of observables.
 * The system dependent variables are stored in World and
 * Config. */

#ifndef __SIMULATION_H_INCLUDED__
#define __SIMULATION_H_INCLUDED__
#include <math.h>
#include <cmath>
#include <vector>
#include <array>
#include <iostream>
#include <string>
#include <typeinfo>
#include <algorithm>
#include <ctime>
#include "World.h"
#include "Config.h"
#include "Particle.h"
#include "Random.h"
#include "Geometry.h"
#include "ParticleType.h"
#include "ParticleInteraction.h"
#include "Reaction.h"
#include "ReactionEvent.h"
#include "UnimolecularCandidate.h"
#include "Neighborlist.h"
#include "utils.h"

#include "Observable.h"
#include "Trajectory.h"
#include "RadialDistribution.h"
#include "MeanSquaredDisplacement.h"
#include "ProbabilityDensity.h"
#include "Energy.h"
#include "Acceptance.h"
#include "ParticleNumbers.h"

class Simulation
{
public:
	Simulation();
	~Simulation();
	/* Start the simulation. Iterate for maxTime timesteps.*/
	void run(const unsigned long maxTime);

	/* World stores positions of particles and other variables
	 * that change during the simulation. Config stores information
	 * that does not change during the simulation. Together they
	 * contain all the system information. */
	World * world;
	Config * config;
	/* Determine if reversibility and/or neighborlist will be used. */
	bool isReversibleDynamics;
	bool isReversibleReactions;
	bool useNeighborList;

	/* Observable methods */
	void writeAllObservablesToFile();
	void writeLastObservableToFile();
	std::string showObservables();
	void deleteAllObservables();
	void deleteLastObservable();
	void new_Trajectory(unsigned long recPeriod, std::string filename);
	void new_RadialDistribution(
		unsigned long recPeriod,
		std::string filename,
		std::vector<double> ranges,
		std::vector< std::vector<unsigned> > considered);
	void new_MeanSquaredDisplacement(
		unsigned long recPeriod,
		std::string filename,
		unsigned particleTypeId);
	void new_ProbabilityDensity(
		unsigned long recPeriod,
		std::string filename,
		unsigned pTypeId,
		std::vector<double> range,
		unsigned int coord);
	void new_Energy(unsigned long recPeriod, std::string filename);
	void new_Acceptance(
		unsigned long recPeriod,
		std::string filename,
		bool reactionsOrDynamics);
	void new_ParticleNumbers(
		unsigned long recPeriod,
		std::string filename,
		unsigned particleTypeId);

private:
	Random * random; // the random number generator
	Utils * utils; // functions for calculating distance
	Neighborlist * neighborlist;
	std::vector< std::unique_ptr<Observable> > observables;

	/*------- core functions/variables not accessible to python -------*/

	/* Configure all observables, e.g. MSD gets all initial coordinates */
	void configureObservables();
	/* Store the objects describing the current state:
	 * energy, particles (positions, forces) and activePairs */
	void saveOldState();
	/* Doing the opposite of the above */
	void restoreOldState();
	/* Perform Brownian Dynamics step on particles
	 * according to their accumulated forces */
	void propagateDynamics();
	/* Perform reactions. Unimolecular according to their
	 * reaction rates and bimolecular only when the pair
	 * is part of activePairs. Already return the ratios
	 * of forward and backward probabilities. */
	double propagateReactions();
	/* Call every observables' record() function, if the timeIndex
	 * corresponds to the predefined recording interval.*/
	void recordObservables(unsigned long int timeIndex);
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
	void calculateSingleForceEnergyCheckReactionCandidate(
		unsigned int indexI,
		unsigned int indexJ);
	void calculateGeometryForcesEnergies();
	void resetForces();
	void resetReactionCandidates();
	double acceptanceDynamics();
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

#endif // __SIMULATION_H_INCLUDED__