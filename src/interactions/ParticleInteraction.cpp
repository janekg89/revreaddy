/* ParticleInteraction.cpp */

#include "ParticleInteraction.h"

ParticleInteraction::ParticleInteraction() { }

ParticleInteraction::~ParticleInteraction() { }

void ParticleInteraction::calculateForceEnergy(
	std::vector<double>& forceI,//out
	double& energy,
	const std::vector<double>& r_ij,//in
	const double& rSquared,
	const double& radiiSquared)
{ }

double ParticleInteraction::calculateEnergy(
	const double rSquared, // in
	const double radiiSquared) // in
{ return 0.; }

/* assuming that affectedTuple is sorted and has length 2*/
bool ParticleInteraction::isAffected(unsigned int i, unsigned int j) {
	if (i <= j) {
		return ( ( i == affectedTuple[0] ) && ( j == affectedTuple[1] ) );
	}
	else {
		return ( ( i == affectedTuple[1] ) && ( j == affectedTuple[0] ) );
	}
	throw Exception("Exceptional state reached in 'isAffected'");
}