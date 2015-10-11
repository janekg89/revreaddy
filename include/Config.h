/* Config.h
 * Contains all information that don't change during run.
 * Such as setting new forces, observables, timestep ... */

#ifndef __CONFIG_H_INCLUDED__
#define __CONFIG_H_INCLUDED__
#include <vector>
#include <string>
#include <iostream>
#include <memory>
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
#include "Fusion2.h"
#include "Fusion3.h"

class Config {
public:
	/* random belongs to Simulation, it is only borrowed and must not
	 * be destroyed. This is because Reactions need random numbers
	 * to calculate probabilites. TODO resolve this */
	Config(World * inWorld, Random * inRandom);
	~Config();
	World * world;
	Random * random;

	std::vector<ParticleType> typeDict;
	/* All forces between particles. Since some reactions need to know
	 * the energies that occur between particles, these are shared
	 * pointers */
	std::vector< std::shared_ptr<ParticleInteraction> > possibleInteractions;
	/* All first order potentials. Used to build geometries. Config owns these
	 * so the pointers are unique */ 
	std::vector< std::unique_ptr<Geometry> > geometries;
	/* Storage of all active reversible reactions that may happen */
	std::vector< std::unique_ptr<Reaction> > possibleReactions;
	/* Stores children of Observable. Those are called in
	 * predefined intervals during run(). */
	std::vector< std::unique_ptr<Observable> > observables;

	unsigned long int maxTime; // length of the simulation
	double timestep;           // the timestep: usually 0.001 or smaller
	double temperature;
	double kBoltzmann;
	bool isPeriodic;           // use periodic boundary conditions or not
	double boxsize;            // length of the periodic simulationbox
	bool isReversibleDynamics;
	bool isReversibleReactions;
	bool useNeighborList;
	unsigned numberBoxes;
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
	std::string getDictName(unsigned i);
	double getDictRadius(unsigned i);
	double getDictDiffusionConstant(unsigned i);
	double getDictReactionRadius(unsigned i);
	unsigned int  getParticleNumber();
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
	std::string getForceName(unsigned i);
	std::string getForceType(unsigned i);
	std::vector<unsigned int> getForceAffectedTuple(unsigned i);
	std::vector<double> getForceParameters(unsigned i);
	double getForceCutoff(unsigned i);
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
	void configure_Fusion(
		unsigned reactionIndex,
		std::vector<unsigned> interactionsIndices,
		double inversePartition,
		double maxDistr,
		double radiiSum,
		double reactionRadiiSum,
		double meanDistr,
		double inverseTemperature);
	void new_Fusion2(
		std::string name,
		unsigned forwardTypeA,
		unsigned forwardTypeB,
		unsigned backwardTypeC,
		double forwardRate,
		double backwardRate);
	void configure_Fusion2(
		unsigned reactionIndex,
		std::vector<unsigned> interactionsIndices,
		double inversePartition,
		double maxDistr,
		double radiiSum,
		double reactionRadiiSum,
		double meanDistr,
		double inverseTemperature);
	void new_Fusion3(
		std::string name,
		unsigned forwardTypeA,
		unsigned forwardTypeB,
		unsigned backwardTypeC,
		double forwardRate,
		double backwardRate);
	void configure_Fusion3(
		unsigned reactionIndex,
		std::vector<unsigned> interactionsIndices,
		double inversePartition,
		double maxDistr,
		double radiiSum,
		double reactionRadiiSum,
		double meanDistr,
		double inverseTemperature,
		double radiusA,
		double radiusB);
	unsigned getNumberReactions();
	std::string getReactionName(unsigned i);
	std::string getReactionType(unsigned i);
	std::vector<unsigned> getReactionForwardTypes(unsigned i);
	std::vector<unsigned> getReactionBackwardTypes(unsigned i);
	double getReactionForwardRate(unsigned i);
	double getReactionBackwardRate(unsigned i);
};

#endif //__CONFIG_H_INCLUDED__