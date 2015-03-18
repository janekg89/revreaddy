/* Force.cpp
 * author: Christoph Froehner
 */

#include "Force.h"
// 2^(1/6), 2^(-1/6) and 2^(-1/3)
#define TWO_POW_ONE_SIXTH 1.122462048309373 
#define TWO_POW_MIN_ONE_SIXTH 0.8908987181403393
#define TWO_POW_MIN_ONE_THIRD 0.7937005259840998

void Force::repulsionForce(
	std::vector<double>& forceI,
	std::vector<double>& r_ij,
	double& rSquared,
	double& radiiSquared,
	double& strength,
	unsigned int& typeI,
	unsigned int& typeJ)
{
	if ( (typeI == 1) && (typeJ == 1) ) {
		double correctedSigmaSquared = TWO_POW_MIN_ONE_THIRD * radiiSquared;
		this->LJ1206Force(
			forceI,
			r_ij,
			rSquared,
			correctedSigmaSquared,
			strength);
		return;
	}
	else {
		this->softcoreForce(forceI, r_ij, rSquared, radiiSquared, strength);
		return;
	}
}

void Force::softcoreForce(
	std::vector<double>& forceI,
	std::vector<double>& r_ij,
	double& rSquared, 
	double& radiiSquared,
	double& strength)
{
	if ( rSquared > radiiSquared ) {
		forceI = {0.,0.,0.};
		return;
	}
	double preFactor = 2. * strength * (1. - sqrt(radiiSquared/rSquared));
	forceI[0] = preFactor * r_ij[0];
	forceI[1] = preFactor * r_ij[1]; 
	forceI[2] = preFactor * r_ij[2];
}

void Force::LJ1206Force(
	std::vector<double>& forceI,
	std::vector<double>& r_ij,
	double& rSquared,
	double& sigmaSquared,
	double& strength)
{
	if ( rSquared > (6.25*sigmaSquared) ) {
		forceI = {0.,0.,0.};
		return;
	}
	double preFactor = -4. * strength;
	preFactor *= pow(sigmaSquared/rSquared,6.) * (12./rSquared) 
			   - pow(sigmaSquared/rSquared,3.) * (6. /rSquared);
	forceI[0] = preFactor * r_ij[0];
	forceI[1] = preFactor * r_ij[1];
	forceI[2] = preFactor * r_ij[2];
}

void Force::repulsionEnergy(
	double& energy,
	double& rSquared,
	double& radiiSquared,
	double& strength,
	unsigned int& typeI,
	unsigned int& typeJ)
{
	if ( (typeI == 1) && (typeJ == 1) ) {
		double correctedSigmaSquared = TWO_POW_MIN_ONE_THIRD * radiiSquared;
		this->LJ1206Energy(
			energy,
			rSquared,
			correctedSigmaSquared,
			strength);
		return;
	}
	else {
		this->softcoreEnergy(
			energy,
			rSquared,
			radiiSquared,
			strength);
		return;
	}
}

void Force::softcoreEnergy(
	double& energy,
	double& rSquared,
	double& radiiSquared,
	double& strength)
{
	if (rSquared > radiiSquared) {
		energy = 0.;
		return;
	}
	// E = strength * ( r - radii )**2
	energy	= sqrt(rSquared) - sqrt(radiiSquared);
	energy *= energy;
	energy *= strength;
}

void Force::LJ1206Energy(
	double& energy,
	double& rSquared,
	double& sigmaSquared,
	double& strength)
{
	if ( rSquared > (6.25*sigmaSquared) ) {
		energy = 0.;
		return;
	}
	energy	= 4. * strength;
	energy *= pow(sigmaSquared/rSquared,6.) - pow(sigmaSquared/rSquared,3.);
}

void Force::repulsionForceEnergy(
	std::vector<double>& forceI,
	double& energy,
	std::vector<double>& r_ij,
	double& rSquared,
	double& radiiSquared,
	double& strength,
	unsigned int& typeI,
	unsigned int& typeJ)
{
	if ( (typeI == 1) && (typeJ == 1) ) {
		double correctedSigmaSquared = TWO_POW_MIN_ONE_THIRD * radiiSquared;
		this->LJ1206ForceEnergy(
			forceI,
			energy,
			r_ij,
			rSquared,
			correctedSigmaSquared,
			strength);
		return;
	}
	else {
		this->softcoreForceEnergy(
			forceI,
			energy,
			r_ij,
			rSquared,
			radiiSquared,
			strength);
		return;
	}
}

void Force::softcoreForceEnergy(
	std::vector<double>& forceI,//out
	double& energy,
	std::vector<double>& r_ij,//in
	double& rSquared,
	double& radiiSquared,
	double& strength)
{
	if ( rSquared > radiiSquared ) {
		forceI = {0.,0.,0.};
		energy = 0.;
		return;
	}
	// E=strength*(r - radii)**2, F=2.*strength*(r - radii)/r * r_ij(vec)
	double a = sqrt(rSquared);
	double c = a - sqrt(radiiSquared);
	energy = strength * c * c;
	a = 2. * strength * c / a;
	forceI[0] = a * r_ij[0];
	forceI[1] = a * r_ij[1];
	forceI[2] = a * r_ij[2];
}

void Force::LJ1206ForceEnergy(
	std::vector<double>& forceI,//out
	double& energy,
	std::vector<double>& r_ij,//in
	double& rSquared,
	double& sigmaSquared,
	double& strength)	
{
	if ( rSquared > (6.25*sigmaSquared) ) {
		energy = 0.;
		forceI = {0.,0.,0.,};
		return;
	}
	double a = pow(sigmaSquared/rSquared, 3.);
	double b = a * a;
	energy = 4. * strength * ( b - a );
	a = -24. * strength * ( 2. * b - a) / rSquared;
	forceI[0] = a * r_ij[0];
	forceI[1] = a * r_ij[1];
	forceI[2] = a * r_ij[2];
}
