/* ParticleType.cpp */

#include "ParticleType.h"

ParticleType::ParticleType(
	std::string inName,
	double inRadius,
	double inDiffusionConstant,
	double inReactionRadius)
{
	this->name = inName;
	this->radius = inRadius;
	this->diffusionConstant = inDiffusionConstant;
	this->reactionRadius = inReactionRadius;
}
