/* Simulation.h
 * author: Christoph Froehner
 *
 * This is the main class, which uses all other classes.
 * Simulation performs the time loop which propagates
 * the particles.
 */

#ifndef __SIMULATION_H_INCLUDED__
#define __SIMULATION_H_INCLUDED__
class Simulation;
#include <math.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <string>
#include <typeinfo>
#include "Particle.h"
#include "Random.h"
#include "Force.h"
#include "Observable.h"
#include "Trajectory.h"
#include "RadialDistribution.h"
#include "MeanSquaredDisplacement.h"
#include "ProbabilityDensity.h"
#include "Energy.h"
#include "Geometry.h"
#include "Wall.h"
#include "DoubleWellZ.h"
#include "TypeDict.h"
#include "utils.h"

class Simulation
{
	public:
		Random * random;                // the random number generator
		Force * force;                  // the force/energy handler
		TypeDict * typeDict;			// dictionary for radii etc.
		std::vector<Particle> activeParticles;
		/* Stores children of Observable */
		std::vector<Observable*> observables;
		/* For later adaptive timestepping methods */
		//std::vector<Particle> consideredParticles; 
		std::vector<Geometry*> geometries;
		unsigned long int maxTime; // length of the simulation
		double timestep;           // the timestep: usually 0.001
		double cumulativeRuntime;  // keeps track of the advanced time
		double temperature;
		double kBoltzmann;
		bool isPeriodic;           // use periodic boundary conditions or not
		double boxsize;            // length of the periodic simulationbox
		double repulsionStrength;  // force constant for particle repulsion
		double energy;
		double oldEnergy;
		unsigned long int acceptions;
		unsigned long int rejections;
		bool isReversible;
		unsigned long long uniqueIdCounter;

		void run();
		void saveOldState();//oldEnergy=energy, oldPos=pos, oldForce=force
		void propagate();
		void recordObservables(unsigned long int timeIndex);
		/* double loop (i,j) over activeParticles and call 
		 * according Forcetype for particle pair (i,j) */
		void calculateRepulsionForcesEnergies(); 
		void calculateGeometryForcesEnergies();
		void resetForces();
		void acceptOrReject();

		Simulation(bool hasDefaultTypes);
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
			double reactionRadius,
			unsigned int forceType);
		std::vector<std::string> getDictNames();
		std::vector<double> getDictRadii();
		std::vector<double> getDictDiffusionConstants();
		std::vector<double> getDictReactionRadii();
		std::vector<unsigned int> getDictForceTypes();
		int  getParticleNumber();
		void deleteAllParticles();
		void writeAllObservablesToFile();
		std::string showObservables();
		void deleteAllObservables();
		void new_Trajectory(std::string filename);
		void new_RadialDistribution(
			std::string filename,
			std::vector<double> ranges,
			std::vector< std::vector<unsigned int> > considered);
		void new_MeanSquaredDisplacement(
			std::string filename,
			unsigned int particleTypeId);
		void new_ProbabilityDensity(
			std::string filename,
			unsigned int pTypeId,
			std::vector<double> range,
			unsigned int coord);
		void new_Energy(unsigned long int recPeriod, std::string filename);
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
};

#endif // __SIMULATION_H_INCLUDED__
