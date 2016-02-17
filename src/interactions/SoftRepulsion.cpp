/* SoftRepulsion.cpp */

#include "SoftRepulsion.h"

// cutoff is set by Simulation, here it is radius_i + radius_j
SoftRepulsion::SoftRepulsion(std::string inName, std::vector <unsigned int> inAffectedTuple, double inRepulsionStrength) {
	this->name = inName;
	this->type = "SoftRepulsion";
	this->parameters = { inRepulsionStrength };
	// apply the convention that the tuple must be sorted
	if ( inAffectedTuple[0] <= inAffectedTuple[1] ) {
		this->affectedTuple.push_back(inAffectedTuple[0]);
		this->affectedTuple.push_back(inAffectedTuple[1]);
	}
	else {
		this->affectedTuple.push_back(inAffectedTuple[1]);
		this->affectedTuple.push_back(inAffectedTuple[0]);	
		LOG_INFO("SoftRepulsion affectedTuple order was inverted")
	}
	this->c = 0.;
	this->a = 0.;
}

void SoftRepulsion::calculateForceEnergy(
	std::vector<double>& forceI, //out
	double& energy,
	const std::vector<double>& r_ij, //in
	const double& rSquared,
	const double& radiiSquared)
{
	if ( rSquared > radiiSquared ) {
		forceI[0]=0.; forceI[1]=0.; forceI[2]=0.;
		energy = 0.;
		return;
	}
	// E=strength*(r - radii)**2, F=2.*strength*(r - radii)/r * r_ij(vec)
	this->a = sqrt(rSquared);
	this->c = this->a - sqrt(radiiSquared);
	energy = this->parameters[0] * this->c * this->c;
	this->a = 2. * this->parameters[0] * this->c / this->a;
	forceI[0] = this->a * r_ij[0];
	forceI[1] = this->a * r_ij[1];
	forceI[2] = this->a * r_ij[2];
}

double SoftRepulsion::calculateEnergy(
	const double rSquared, //in
	const double radiiSquared) //in
{
	if ( rSquared > radiiSquared ) {
		return 0.;
	}
	// E = strength * (r - radii)**2
	this->c = sqrt(rSquared);
	this->c -= sqrt(radiiSquared);
	return this->parameters[0] * this->c * this->c;
}