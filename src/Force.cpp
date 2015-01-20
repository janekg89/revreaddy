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
	std::string& typeI,
	std::string& typeJ)
{
	if ( (typeI == "soft") && (typeJ == "soft") ) {
		this->softcoreForce(forceI, r_ij, rSquared, radiiSquared, strength);
		return;
	}
	if ( (typeI == "lj") && (typeJ == "lj") ) {
		double correctedSigmaSquared = TWO_POW_MIN_ONE_THIRD * radiiSquared;
		this->LJ1206Force(forceI, r_ij, rSquared, correctedSigmaSquared, strength);
		return;
	}
	if ( (typeI == "lj") && (typeJ == "soft") ) {
		double correctedSigmaSquared = TWO_POW_MIN_ONE_THIRD * radiiSquared;
		this->LJ1206Force(forceI, r_ij, rSquared, correctedSigmaSquared, strength);
		return;
	}
	if ( (typeJ == "soft") && (typeJ == "lj") ) {
		double correctedSigmaSquared = TWO_POW_MIN_ONE_THIRD * radiiSquared;
		this->LJ1206Force(forceI, r_ij, rSquared, correctedSigmaSquared, strength);
		return;
	}
	else {
		throw "Particle types are neither 'soft' nor 'lj' !";
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
	double preFactor = strength * (1. - sqrt(radiiSquared/rSquared));
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
