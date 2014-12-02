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
#include "Particle.h"
#include "Random.h"
#include "Potential.h"
#include "Observable.h"
#include <array>
#include <vector>
#include <iostream>

class Simulation
{
	public:
		Random * random; 				// the random number generator
		Potential * potential; 			// the force/energy handler
		std::vector<Particle> activeParticles;
		std::vector<Observable*> observables;	// stores children of Observable
		//std::vector<Particle> consideredParticles; // for later adaptive timestepping methods
		unsigned long int maxTime;		// length of the simulation
		double timestep;				// the timestep: usually 0.001
		double temperature;
		double kBoltzmann;
		bool isPeriodic;				// use periodic boundary conditions or not
		double boxsize;					// length of the periodic simulationbox
		double repulsionStrength; 		// force constant for softcore particle repulsion

		void addParticle(std::array<double,3> initPos, double rad, double diffConst);
		void run();
		void propagate();
		void recordObservables(unsigned long int t);
		// should double loop (i,j) over activeParticles and call according Forcetype for
		// particle pair (i,j)
		void calculateRepulsionForces();
		std::array<double,3> getMinDistance(std::array<double,3> r_i, std::array<double,3> r_j);

		Simulation();
		virtual ~Simulation();

};

#endif
