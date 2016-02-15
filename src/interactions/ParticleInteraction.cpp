/* ParticleInteraction.cpp */

#include "ParticleInteraction.h"
#define print(a) std::cout << a << std::endl;

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

double ParticleInteraction::calculateEnergy(
	double rSquared, // in
	double radiiSquared) // in
{
	return 0.;
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
	print("Error: exceptional state reached in 'isAffected'")
	return false;
}