/* ParticleType.h
 * An object describing the properties of one particle type.
 * The particleType id is determined by its position in the
 * vector typeDict in the main class Simulation.*/

#ifndef __PARTICLETYPE_H_INCLUDED__
#define __PARTICLETYPE_H_INCLUDED__
#include <vector>
#include <string>

class ParticleType
{
	public:
		std::string name;
		double radius;
		double diffusionConstant;
		double reactionRadius;

		ParticleType(
			std::string inName,
			double inRadius,
			double inDiffusionConstant,
			double inReactionRadius);
};

#endif //__PARTICLETYPE_H_INCLUDED__
