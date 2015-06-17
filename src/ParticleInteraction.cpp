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

/* assuming that affectedTuple is sorted and has length 2*/
bool ParticleInteraction::isAffected(unsigned int i, unsigned int j)
{
	if (i <= j) {
		return ( ( i == affectedTuple[0] ) && ( j == affectedTuple[1] ) );
	}
	else {
		return ( ( i == affectedTuple[1] ) && ( j == affectedTuple[0] ) );
	}
	std::cout << "Error: exceptional state reached in 'isAffected'\n";
	return false;
}
