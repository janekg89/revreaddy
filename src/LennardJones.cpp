/* LennardJones.cpp */

#include "LennardJones.h"
// 2^(-1/3)
#define TWO_POW_MIN_ONE_THIRD 0.7937005259840998

//cutoff is set by Simulation, here it is 2.5*(radius_i+radius_j)
LennardJones::LennardJones(
	std::string inName,
	std::vector<unsigned int> inAffectedTuple,
	double inEpsilon)
{
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
		std::cout << "Info: LennardJones affectedTuple order was inverted\n";
	}
}

void LennardJones::calculateForceEnergy(
	std::vector<double>& forceI, //out
	double& energy,
	std::vector<double>& r_ij, //in
	double& rSquared,
	double& radiiSquared)
{
	double sigmaSquared = TWO_POW_MIN_ONE_THIRD * radiiSquared;
	if ( rSquared > (6.25*sigmaSquared) ) {
		energy = 0.;
		forceI = {0.,0.,0.,};
		return;
	}
	double a = pow(sigmaSquared/rSquared, 3.);
	double b = a * a;
	energy = 4. * this->parameters[0] * ( b - a );
	a = -24. * this->parameters[0] * ( 2. * b - a) / rSquared;
	forceI[0] = a * r_ij[0];
	forceI[1] = a * r_ij[1];
	forceI[2] = a * r_ij[2];
}
