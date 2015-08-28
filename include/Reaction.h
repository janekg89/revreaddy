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
#include "Simulation.h"

class Reaction
{
	public:
		std::string name;
		std::vector<unsigned int> forwardTypes;
		std::vector<unsigned int> backwardTypes;
		double forwardRate;
		double backwardRate;
		
		/* The following perform the reaction in the
		 * given direction. particleIndices contains the 
		 * activeParticles index of particles on which the reaction
		 * will be performed. This is either one or 
		 * two numbers, since reactions are either 
		 * uni or bimolecular. */
		/* NOTE: These functions manipulate Simulation.activeParticles
		 * but should do this only by using the methods provided by
		 * Simulation, e.g. addParticle() */
		virtual double performForward(
			std::vector<unsigned long int> particleIndices);
		virtual double performBackward(
			std::vector<unsigned long int> particleIndices);
};

#endif //__REACTION_H_INCLUDED__
