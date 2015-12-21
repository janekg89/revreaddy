/* Simulation.h
 * author: Christoph Froehner
 *
 * This is the main class, which uses all other classes.
 * Simulation performs the time loop which propagates
 * the particles. */

#ifndef __SIMULATION_H_INCLUDED__
#define __SIMULATION_H_INCLUDED__
/* This forward declaration is still necessary due 
 * to a non-resolved circular dependence with one or more
 * observables and Reactions, since they manipulate the
 * current state directly. */
class Simulation;
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
#include "Observable.h"
#include "Geometry.h"
#include "ParticleType.h"
#include "ParticleInteraction.h"
#include "Reaction.h"
#include "ReactionEvent.h"
#include "UnimolecularCandidate.h"
#include "Neighborlist.h"
#include "utils.h"

class Simulation
{
	public:
		Simulation();
		~Simulation();

		World * world;
		Config * config;
		Random * random; // the random number generator
		Utils * utils;

		Neighborlist * neighborlist;

		/*------- core functions/variables not accessible to python -------*/

		/* Store the objects describing the current state:
		 * energy, activeParticles (positions, forces) and activePairs */
		void saveOldState();
		/* Doing the opposite of the above */
		void restoreOldState();
		/* Perform Brownian Dynamics step on activeParticles
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
		/* double loop (i,j) over activeParticles and call 
		 * according Forcetype for particle pair (i,j) --> O(n^2) */
		void calculateInteractionForcesEnergiesNaive();
		/* does the same as above, but only considers interactions
		 * of neighboring boxes, that have the size of the maximum
		 * cutoff distance --> O(n) */
		void calculateInteractionForcesEnergiesWithLattice(unsigned int numberBoxes);
		/* evaluate the force and energy for a given pair of
		 * particles and store their unique ids in activePairs
		 * if they are in reactive distance */
		void calculateSingleForceEnergy(
			unsigned int indexI,
			unsigned int indexJ);
		void calculateSingleForceEnergyOnlyForI(
			unsigned int indexI,
			unsigned int indexJ);
		void calculateGeometryForcesEnergies();
		void resetForces();
		void resetActivePairs();
		double acceptanceDynamics();
		double acceptanceReactions();
		bool acceptOrReject(double acceptance);
		/* Return the position in activeParticles of the 
		 * particle with the uniqueId id. Particles are added
		 * to activeParticles only using addParticle(), which
		 * means that it is sorted w.r.t. uniqueIds. Therefore
		 * this function performs a binary search with
		 * complexity O(log n). The return value is signed so
		 * the case "no particle found" is expressed by "-1" */
		long findParticleIndex(unsigned long long id);

		/* Start the simulation. Iterate for maxTime timesteps.*/
		void run();

		/* members that are once allocated so that less allocation
		 * appears on the run */
		std::vector<double> forceI;
		std::vector<double> forceJ;
		double energyBuffer;
		std::vector<double> r_ij;	
		double rSquared;
		double radiusI;
		double radiusJ;
		double radiiSquared;
		double reactionRadiiSquared;
		unsigned sizePossibleInteractions;
};

#endif // __SIMULATION_H_INCLUDED__