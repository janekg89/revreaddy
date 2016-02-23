/* ParticleType.cpp */

#include "ParticleType.h"

ParticleType::ParticleType(
	std::string inName,
	double inRadius,
	double inDiffusionConstant)
{
	this->name = inName;
	this->radius = inRadius;
	this->diffusionConstant = inDiffusionConstant;
}
