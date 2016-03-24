/* LennardJones.cpp */

#include "LennardJones.h"
// 2^(-1/3)
#define TWO_POW_MIN_ONE_THIRD 0.7937005259840998

//cutoff is set by Simulation, here it is 2.5*(radius_i+radius_j)
LennardJones::LennardJones(std::string inName, std::array<unsigned,2> inAffectedTuple, double inEpsilon) {
	this->name = inName;
	this->type = "LennardJones";
	this->parameters = { inEpsilon };
	this->affectedTuple = inAffectedTuple;
	// apply the convention that the tuple must be sorted
	std::sort( affectedTuple.begin(), affectedTuple.end() );
}

void LennardJones::calculateForceEnergy(
	std::vector<double>& forceI, //out
	double& energy,
	const std::vector<double>& r_ij, //in
	const double& rSquared,
	const double& radiiSquared)
{
	// radiiSquared = (R_i + R_j)^2, R_i - particle radius of i
	// rSquared = r_ij^2
	if ( rSquared > (6.25*radiiSquared) ) {
		energy = 0.;
		forceI = {0.,0.,0.,};
		return;
	}
	double sigmaSquared = TWO_POW_MIN_ONE_THIRD * radiiSquared;
	double a = pow(sigmaSquared/rSquared, 3.);
	double b = a * a;
	energy = 4. * this->parameters[0] * ( b - a );
	a = -24. * this->parameters[0] * ( 2. * b - a) / rSquared;
	forceI[0] = a * r_ij[0];
	forceI[1] = a * r_ij[1];
	forceI[2] = a * r_ij[2];
}

double LennardJones::calculateEnergy(
	const double rSquared,
	const double radiiSquared)
{
	if ( rSquared > (6.25*radiiSquared) ) {
		return 0.;
	} 
	double sigmaSquared = TWO_POW_MIN_ONE_THIRD * radiiSquared;
	double a = pow(sigmaSquared/rSquared, 3.);
	double b = a * a;
	return 4. * this->parameters[0] * ( b - a );
}