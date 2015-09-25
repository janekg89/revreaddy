/* Reaction.h
 * The Reaction object is instantiated for every type of
 * reaction. More precisely, this class will function as
 * the parent type and actual reaction types (conversion,
 * decay, birth, enzymatic) are implemented by children
 * of this class. */

#ifndef __REACTION_H_INCLUDED__
#define __REACTION_H_INCLUDED__
#include <string>
#include <vector>
#include "World.h"
#include "Particle.h" // positions must be obtained

class Reaction
{
	public:
		Reaction();
		virtual ~Reaction();
		std::string name;
		std::vector<unsigned int> forwardTypes;
		std::vector<unsigned int> backwardTypes;
		double forwardRate;
		double backwardRate;
		
		/* Check if the given types or type are/is 
		 * affected by the reaction in forward direction. */
		bool isAffectedForward(std::vector<unsigned int> types);
		bool isAffectedForward(unsigned int type);
		/* Same as above but checks for backward direction */
		bool isAffectedBackward(std::vector<unsigned int> types);
		bool isAffectedBackward(unsigned int type);
		/* The following perform the reaction in the
		 * given direction. particleIndices contains the 
		 * activeParticles index of particles on which the reaction
		 * will be performed. This is either one or 
		 * two numbers, since reactions are either 
		 * uni or bimolecular. */
		/* NOTE: These functions manipulate world->activeParticles
		 * but should do this only by using the methods provided by
		 * World, e.g. add/removeParticle(). That is why they receive
		 * a pointer to the Simulation that called it. */
		/* NOTE2: particleIndices are the positions in activeParticles 
		 * not the uniqueIds. */
		virtual double performForward(
			std::vector<unsigned long int> particleIndices,
			World * world);
		virtual double performBackward(
			std::vector<unsigned long int> particleIndices,
			World * world);
};

#endif //__REACTION_H_INCLUDED__