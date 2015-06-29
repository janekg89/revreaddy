/* Simulation.h
 * author: Christoph Froehner
 *
 * This is the main class, which uses all other classes.
 * Simulation performs the time loop which propagates
 * the particles.
 */

#ifndef __SIMULATION_H_INCLUDED__
#define __SIMULATION_H_INCLUDED__
/* This forward declaration is still necessary due 
 * to a non-resolved circular dependence with one or more
 * observables. */
class Simulation;
#include <math.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <string>
#include <typeinfo>
#include <algorithm>
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
#include "TypeDict.h"
#include "ParticleInteraction.h"
#include "SoftRepulsion.h"
#include "LennardJones.h"
#include "utils.h"

class Simulation
{
	public:
		Random * random;                // the random number generator
		TypeDict * typeDict;			// dictionary for radii etc.
		std::vector<Particle> activeParticles;
		/* Stores children of Observable */
		std::vector<Observable*> observables;
		/* For later adaptive timestepping methods */
		//std::vector<Particle> consideredParticles; 
		std::vector<Geometry*> geometries;
		/* Storage of all active reactions that may happen */
		//std::vector<Reaction> reactions;
		/* All forces between particles */
		std::vector<ParticleInteraction*> possibleForces;

		unsigned long int maxTime; // length of the simulation
		double timestep;           // the timestep: usually 0.001 or smaller
		double cumulativeRuntime;  // keeps track of the advanced time
		double temperature;
		double kBoltzmann;
		bool isPeriodic;           // use periodic boundary conditions or not
		double boxsize;            // length of the periodic simulationbox
		double energy;
		double oldEnergy;
		double currentAcceptance;  // the last calculated acceptance prob
		unsigned long int acceptions;
		unsigned long int rejections;
		bool isReversible;
		unsigned long long uniqueIdCounter;
		bool useNeighborList;

		void run();
		void saveOldState();//oldEnergy=energy, oldPos=pos, oldForce=force
		void propagate();
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
		/* evaluate the force and energy for
		 * a given pair of particles */
		void calculateSingleForceEnergy(
			unsigned int indexI,
			unsigned int indexJ);
		void calculateGeometryForcesEnergies();
		void resetForces();
		void acceptOrReject();

		Simulation();
		~Simulation();

		/*------- functions that will be wrapped by python -----------*/

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
		std::vector<std::string> getDictNames();
		std::vector<double> getDictRadii();
		std::vector<double> getDictDiffusionConstants();
		std::vector<double> getDictReactionRadii();
		int  getParticleNumber();
		void deleteAllParticles();
		void writeAllObservablesToFile();
		std::string showObservables();
		void deleteAllObservables();
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
};

#endif // __SIMULATION_H_INCLUDED__
