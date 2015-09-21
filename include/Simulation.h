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
#include <iostream>
#include <string>
#include <typeinfo>
#include <algorithm>
#include <time.h>
#include "World.h"
#include "Config.h"
#include "Particle.h"
#include "Random.h"
#include "Observable.h"
#include "Trajectory.h"
#include "RadialDistribution.h"
#include "MeanSquaredDisplacement.h"
#include "ProbabilityDensity.h"
#include "Energy.h"
#include "Acceptance.h"
#include "Geometry.h"
#include "Wall.h"
#include "DoubleWellZ.h"
#include "ParticleType.h"
#include "ParticleInteraction.h"
#include "SoftRepulsion.h"
#include "LennardJones.h"
#include "Reaction.h"
#include "Conversion.h"
#include "ReactionEvent.h"
#include "utils.h"

class Simulation
{
	public:
		World * world;
		Config * config;
		Random * random;                // the random number generator

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
		long int findParticleIndex(unsigned long long id);

		/*------- functions that will be wrapped by python -------*/

		Simulation();
		~Simulation();
		/* Start the simulation. Iterate for maxTime timesteps.*/
		void run();
		void addParticle(
			std::vector<double> initPos,
			unsigned int particleTypeId);
		/* Obtain the position of the particle activeParticles[index] */
		std::vector<double> getPosition(int index);
		void                setPosition(int index, std::vector<double> newPos);
		unsigned int getTypeId(int index);
		void         setTypeId(int index, unsigned int typeId);
		void deleteAllParticles();
};

#endif // __SIMULATION_H_INCLUDED__
