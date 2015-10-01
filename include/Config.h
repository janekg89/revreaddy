/* Config.h
 * Contains all information that don't change during run.
 * Such as setting new forces, observables, timestep ... */

#ifndef __CONFIG_H_INCLUDED__
#define __CONFIG_H_INCLUDED__
#include <vector>
#include <string>
#include <iostream>
#include <typeinfo>
#include "World.h"
#include "Random.h"
#include "Observable.h"
#include "Trajectory.h"
#include "RadialDistribution.h"
#include "MeanSquaredDisplacement.h"
#include "ProbabilityDensity.h"
#include "Energy.h"
#include "Acceptance.h"
#include "ParticleNumbers.h"
#include "Geometry.h"
#include "Wall.h"
#include "DoubleWellZ.h"
#include "ParticleType.h"
#include "ParticleInteraction.h"
#include "SoftRepulsion.h"
#include "LennardJones.h"
#include "Reaction.h"
#include "Conversion.h"
#include "Fusion.h"

class Config
{
	public:
		/* random belongs to Simulation, it is only borrowed and must not
		 * be destroyed. This is because Reactions need random numbers
		 * to calculate probabilites. TODO resolve this */
		Config(World * inWorld, Random * inRandom);
		~Config();
		World * world;
		Random * random;

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

		/*------------------------------------------------------------------*/
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
		void new_ParticleNumbers(
			unsigned long recPeriod,
			std::string filename,
			unsigned particleTypeId);
		void deleteAllGeometries();
		void new_Wall(
			std::vector<double> normal,
			std::vector<double> point,
			double strength,
			std::vector<unsigned int> particleTypeIds);
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
		void deleteAllReactions();
		void new_Conversion(
			std::string name,
			unsigned forwardType,
			unsigned backwardType,
			double forwardRate,
			double backwardRate);	
		void new_Fusion(
			std::string name,
			unsigned forwardTypeA,
			unsigned forwardTypeB,
			unsigned backwardTypeC,
			double forwardRate,
			double backwardRate);
};

#endif //__CONFIG_H_INCLUDED__