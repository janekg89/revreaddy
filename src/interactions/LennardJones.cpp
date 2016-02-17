/* LennardJones.cpp */

#include "LennardJones.h"
// 2^(-1/3)
#define TWO_POW_MIN_ONE_THIRD 0.7937005259840998

//cutoff is set by Simulation, here it is 2.5*(radius_i+radius_j)
LennardJones::LennardJones(std::string inName, std::vector<unsigned int> inAffectedTuple, double inEpsilon) {
	this->name = inName;
	this->type = "LennardJones";
	this->parameters = { inEpsilon };
	// apply the convention that the tuple must be sorted
	if ( inAffectedTuple[0] <= inAffectedTuple[1] ) {
		this->affectedTuple.push_back(inAffectedTuple[0]);
		this->affectedTuple.push_back(inAffectedTuple[1]);
	}
	else {
		this->affectedTuple.push_back(inAffectedTuple[1]);
		this->affectedTuple.push_back(inAffectedTuple[0]);	
		LOG_INFO("LennardJones affectedTuple order was inverted")
	}
}

void LennardJones::calculateForceEnergy(
	std::vector<double>& forceI, //out
	double& energy,
	const std::vector<double>& r_ij, //in
	const double& rSquared,
	const double& radiiSquared)
{
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