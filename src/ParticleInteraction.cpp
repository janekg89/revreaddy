/* ParticleInteraction.cpp */

#include "ParticleInteraction.h"

ParticleInteraction::ParticleInteraction()
{

}

ParticleInteraction::~ParticleInteraction()
{

}

void ParticleInteraction::calculateForceEnergy(
	std::vector<double>& forceI,//out
	double& energy,
	std::vector<double>& r_ij,//in
	double& rSquared,
	double& radiiSquared)
{

}

bool ParticleInteraction::isAffected(unsigned int i, unsigned int j)
{
	return false;
}
