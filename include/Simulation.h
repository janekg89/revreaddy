/* Simulation.h
 * author: Christoph Froehner
 *
 * This is the main class, which uses all other classes.
 * Simulation performs the time loop which propagates
 * the particles.
 */

#ifndef __SIMULATION_H_INCLUDED__
#define __SIMULATION_H_INCLUDED__
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
#include "TrajectorySingle.h"
#include "RadialDistribution.h"
#include "MeanSquaredDisplacement.h"
#include "utils.h"

class Simulation
{
	public:
		Random * random;                // the random number generator
		Force * force;                  // the force/energy handler
		std::vector<Particle> activeParticles;
		/* Stores children of Observable */
		std::vector<Observable*> observables;
		/* For later adaptive timestepping methods */
		//std::vector<Particle> consideredParticles; 
		unsigned long int maxTime; // length of the simulation
		double timestep;           // the timestep: usually 0.001
		double cumulativeRuntime;  // keeps track of the advanced time
		double temperature;
		double kBoltzmann;
		bool isPeriodic;           // use periodic boundary conditions or not
		double boxsize;            // length of the periodic simulationbox
		double repulsionStrength;  // force constant for particle repulsion
		bool verbose;
		double energy;
		double oldEnergy;
		unsigned long int acceptions;
		unsigned long int rejections;
		bool isReversible;
		unsigned long long uniqueIdCounter;

		void run();
		void saveOldState();//oldEnergy=energy, oldPos=pos, oldForce=force
		void propagate();
		void recordObservables(double t);
		/* double loop (i,j) over activeParticles and call 
		 * according Forcetype for particle pair (i,j) */
		void calculateRepulsionForcesEnergies(); 
		void resetForces();
		void acceptOrReject();

		Simulation();
		~Simulation();

		/*------- functions that will be wrapped by python -----------*/

		void addParticle(
			std::vector<double> initPos,
			unsigned int particleTypeId,
			double rad,
			double diffConst);

		/* Obtain the position of the particle activeParticles[index] */
		std::vector<double> getPosition(int index);
		void                setPosition(int index, std::vector<double> newPos);
		unsigned int getTypeId(int index);
		void         setTypeId(int index, unsigned int typeId);
		double getRadius(int index);
		void   setRadius(int index, double radius);
		double getDiffusionConstant(int index);
		void   setDiffusionConstant(int index, double diffusionConstant);
		int  getParticleNumber();
		void deleteAllParticles();
		void writeAllObservablesToFile();
		std::string showObservables();
		void deleteAllObservables();
		void new_Trajectory(std::string filename);
		void new_TrajectorySingle();
		std::vector< std::vector<double> > getTrajectorySingle();
		void new_RadialDistribution(
			std::string filename,
			std::vector<double> ranges,
			std::vector< std::vector<unsigned int> > considered);
		void new_MeanSquaredDisplacement(
			std::string filename,
			unsigned int particleTypeId);
};

#endif // __SIMULATION_H_INCLUDED__
