/* World.h
 * This object represents the current state of the simulation
 * and contains mostly only particle positions as a function
 * of time. It is mostly a data structure with modifiers. */

#ifndef __WORLD_H_INCLUDED__
#define __WORLD_H_INCLUDED__
#include <vector>
#include "Particle.h"

class World
{
	public:
		World();
		/* activeParticles has an entry for every particle
		 * which holds its position and typeId.
		 * oldActiveParticles temporarily saves the particle
		 * states in case the timestep is rejected and the old
		 * state must be restored. */
		std::vector<Particle> activeParticles;
		std::vector<Particle> oldActiveParticles;
		/* activePairs has an entry for every pair (i,j) particles'
		 * indices (NOT uniqueIds, since particle indices are conserved
		 * during until reactions are propagated)
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

		void addParticle(
			std::vector<double> initPos,
			unsigned int particleTypeId);
		void removeParticle(unsigned long int index);
		std::vector<double> getPosition(unsigned long int index);
		void setPosition(unsigned long int index, std::vector<double> newPos);
		unsigned int getTypeId(unsigned long int index);
		void setTypeId(unsigned long int index, unsigned int typeId);
		void deleteAllParticles();	
};

#endif //__WORLD_H_INCLUDED__