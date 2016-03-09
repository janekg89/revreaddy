/* Config.h
 * Contains all information regarding the system 
 * that don't change during run.
 * Such as setting new forces, timestep and reactions */

#ifndef __CONFIG_H_INCLUDED__
#define __CONFIG_H_INCLUDED__
#include <vector>
#include <string>
#include <iostream>
#include <memory>
#include <typeinfo>
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
#include "Exception.h"
#include "logging.h"

class Config {
public:
	Config();
	~Config();
	/* The list of particle types, by construction ordered and accessed via
	 * its id (position in the vector) */
	std::vector<ParticleType> particleTypes;
	/* All forces between particles. Since some reactions need to know
	 * the energies that occur between particles, these are shared
	 * pointers */
	std::vector< std::shared_ptr<ParticleInteraction> > interactions;
	/* All first order potentials. Used to build geometries. Config owns these
	 * so the pointers are unique */ 
	std::vector< std::unique_ptr<Geometry> > geometries;
	/* Storage of all active reversible reactions that may happen */
	std::vector< std::unique_ptr<Reaction> > reactions;

	double timestep;
	/* The temperature. When using energy units of kJ/mol then kT = 2.437 at 293 K*/
	double kT;
	bool isPeriodic; // use periodic boundary conditions or not
	double boxsize; // length of the periodic simulationbox

	/*------------------------------------------------------------------*/
	void deleteAllParticleTypes();
	void new_Type(const std::string name, const double radius, const double diffusionConstant);
	unsigned getNumberOfParticleTypes();
	std::string getParticleTypeName(unsigned i);
	double getParticleTypeRadius(unsigned i);
	double getParticleTypeDiffusionConstant(unsigned i);

	void deleteAllGeometries();
	void new_Wall(std::vector<double> normal, std::vector<double> point, double strength, std::vector<unsigned int> particleTypeIds);
	void new_DoubleWellZ(double distanceMinima,	double strength, std::vector<unsigned int> particleTypeIds);

	void deleteAllInteractions();
	void new_SoftRepulsion(std::string name, std::vector<unsigned int> affectedTuple, double repulsionStrength);
	void new_LennardJones(std::string name,	std::vector<unsigned int> affectedTuple, double epsilon);
	unsigned getNumberInteractions();
	std::string getInteractionName(unsigned i);
	std::string getInteractionType(unsigned i);
	std::vector<unsigned int> getInteractionAffectedTuple(unsigned i);
	std::vector<double> getInteractionParameters(unsigned i);
	double getInteractionCutoff(unsigned i);

	void configureReactions();
	void deleteAllReactions();
	void new_Conversion(std::string name, unsigned forwardType, unsigned backwardType, double forwardRate, double backwardRate);	
	void new_Fusion(
		std::string name,
		unsigned forwardTypeA,
		unsigned forwardTypeB,
		unsigned backwardTypeC,
		double forwardRate,
		double backwardRate,
		double reactionDistance);
	void configureFusion(
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
private:
	// TODO
	std::vector<unsigned> sort(std::vector<unsigned> x);
};

#endif //__CONFIG_H_INCLUDED__