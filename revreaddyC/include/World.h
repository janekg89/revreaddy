/* World.h
 * This object represents the current state of the simulation
 * and contains mostly only particle positions as a function
 * of time. It is mostly a data structure with modifiers. */

#ifndef __WORLD_H_INCLUDED__
#define __WORLD_H_INCLUDED__
#include <vector>
#include <boost/multi_array.hpp>
#include <map>
#include "Particle.h"
#include "ReactionEvent.h"
#include "logging.h"
#include "Exception.h"
#include "Random.h"

class World {
public:
	World();
	/* particles has an entry for every particle
	 * which holds its position and typeId.
	 * oldParticles temporarily saves the particle
	 * states in case the timestep is rejected and the old
	 * state must be restored. */
	std::vector<Particle> particles;
	std::vector<Particle> oldParticles;
	/* reactionCandidates holds one ReactionEvent for every
	 * particle or particle-pair that might undergo a reaction.
	 * The particle-pair candidates are found during
	 * calculation of forces. The single particle candidates
	 * are added afterwards. */
	std::vector<ReactionEvent> reactionCandidates;
	std::vector<ReactionEvent> oldReactionCandidates;
	/* energy is the sum of all energy terms in the system by
	 * particle interactions and geometries, i.e. first and second
	 * order potential energies */
	double energy;
	double oldEnergy;
	double cumulativeRuntime; // keeps track of the advanced time
	/* The last calculated acceptance probability for dynamical step */
	double acceptProbDiffusion;
	/* The last calculated acceptance probability for reaction step */
	double acceptProbReactions;
	unsigned long long uniqueIdCounter;
	unsigned long int acceptionsDiffusion;
	unsigned long int rejectionsDiffusion;
	unsigned long int acceptionsReactions;
	unsigned long int rejectionsReactions;

	void addParticle(std::vector<double> initPos, unsigned particleTypeId);
	void removeParticle(unsigned long index);

	unsigned long getNumberOfParticles();
	std::vector<double> getPosition(unsigned long index);
	void setPosition(unsigned long index, std::vector<double> newPos);
	unsigned int getTypeId(unsigned long index);
	void setTypeId(unsigned long index, unsigned typeId);
	unsigned long long getUniqueId(unsigned long index);
	void deleteAllParticles();

	/*Fractional stuff*/
	std::map<unsigned, boost::multi_array<double,2>> increments;
	std::map<unsigned, size_t> incrementsIndex;
	/* Change addParticle and removeParticle, update increments and incrementsIndex*/
	// generate new increments
	void addParticleAndIncrements(std::vector<double> initPos, unsigned particleTypeId, Random * random,
								  unsigned long maxTime, double &timestep, double &diffConst, double &alpha);
	 // erase entries in increments and incrementsIndex
	void removeParticleAndIncrements(unsigned long index);
};

#endif //__WORLD_H_INCLUDED__