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
		/*------- variables for initial setup -------*/
		
		Random * random;                // the random number generator
		std::vector<ParticleType> typeDict;
		/* All forces between particles */
		std::vector<ParticleInteraction*> possibleInteractions;
		/* All first order potentials. Used to build geometries. */ 
		std::vector<Geometry*> geometries;
		/* Storage of all active reversible reactions that may happen */
		std::vector<Reaction*> possibleReactions;
		/* Stores children of Observable. Those are called in
		 * predefined intervals during run(). */
		std::vector<Observable*> observables;

		unsigned long int maxTime; // length of the simulation
		double timestep;           // the timestep: usually 0.001 or smaller
		double temperature;
		double kBoltzmann;
		bool isPeriodic;           // use periodic boundary conditions or not
		double boxsize;            // length of the periodic simulationbox
		bool isReversibleDynamics;
		bool isReversibleReactions;
		bool useNeighborList;
		/* This variable is used only by the author to test different
		 * methods of performing reactions. */
		unsigned int reactionPropagation;

		/*------- state variables, change during run() -------*/

		/* activeParticles has an entry for every particle
		 * which holds its position and typeId.
		 * oldActiveParticles temporarily saves the particle
		 * states in case the timestep is rejected and the old
		 * state must be restored. */
		std::vector<Particle> activeParticles;
		std::vector<Particle> oldActiveParticles;
		/* activePairs has an entry for every pair (i,j) unique ids
		 * that is within reaction range at the current time.
		 * oldActivePairs temporarily saves the state of
		 * activePairs and is restored upon timestep-rejection. */
		std::vector< std::vector<unsigned long> > activePairs;
		std::vector< std::vector<unsigned long> > oldActivePairs;
		/* MAYBE For later adaptive timestepping methods */
		//std::vector<Particle> consideredParticles;
		/* energy is the sum of all energy terms in the system by
		 * particle interactions and geometries, i.e. first and second
		 * order potential energies */
		double energy;
		double oldEnergy;
		double cumulativeRuntime;  // keeps track of the advanced time
		/* The last calculated acceptance probability for dynamical step */
		double acceptProbDynamics;
		/* The last calculated acceptance probability for reaction step */
		double acceptProbReactions;
		unsigned long long uniqueIdCounter;
		unsigned long int acceptionsDynamics;
		unsigned long int rejectionsDynamics;
		unsigned long int acceptionsReactions;
		unsigned long int rejectionsReactions;

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
		void new_Type(
			std::string name,
			double radius,
			double diffusionConstant,
			double reactionRadius);
		unsigned int getNumberOfTypes();
		std::string getDictName(unsigned int i);
		double getDictRadius(unsigned int i);
		double getDictDiffusionConstant(unsigned int i);
		double getDictReactionRadius(unsigned int i);
		unsigned int  getParticleNumber();
		void deleteAllParticles();
		void writeAllObservablesToFile();
		void writeLastObservableToFile();
		std::string showObservables();
		void deleteAllObservables();
		void deleteLastObservable();
		void new_Trajectory(unsigned long int recPeriod, std::string filename);
		void new_RadialDistribution(
			unsigned long int recPeriod,
			std::string filename,
			std::vector<double> ranges,
			std::vector< std::vector<unsigned int> > considered);
		void new_MeanSquaredDisplacement(
			unsigned long int recPeriod,
			std::string filename,
			unsigned int particleTypeId);
		void new_ProbabilityDensity(
			unsigned long int recPeriod,
			std::string filename,
			unsigned int pTypeId,
			std::vector<double> range,
			unsigned int coord);
		void new_Energy(unsigned long int recPeriod, std::string filename);
		void new_Acceptance(unsigned long int recPeriod,std::string filename);
		void deleteAllGeometries();
		void new_Wall(
			std::vector<double> normal,
			std::vector<double> point,
			double strength,
			std::vector<unsigned int>& particleTypeIds);
		void new_DoubleWellZ(
			double distanceMinima,
			double strength,
			std::vector<unsigned int> particleTypeIds);
		void deleteAllForces();
		void new_SoftRepulsion(
			std::string name,
			std::vector<unsigned int> affectedTuple,
			double repulsionStrength);
		void new_LennardJones(
			std::string name,
			std::vector<unsigned int> affectedTuple,
			double epsilon);
		unsigned int getNumberForces();
		std::string getForceName(unsigned int i);
		std::string getForceType(unsigned int i);
		std::vector<unsigned int> getForceAffectedTuple(unsigned int i);
		std::vector<double> getForceParameters(unsigned int i);
		void new_Conversion(
			std::string name,
			std::vector<unsigned int> forwardTypes,
			std::vector<unsigned int> backwardTypes,
			double forwardRate,
			double backwardRate);
};

#endif // __SIMULATION_H_INCLUDED__
