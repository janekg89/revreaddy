/* Force.cpp
 * author: Christoph Froehner
 */

#include "Force.h"
// 2^(1/6), 2^(-1/6) and 2^(-1/3)
#define TWO_POW_ONE_SIXTH 1.122462048309373 
#define TWO_POW_MIN_ONE_SIXTH 0.8908987181403393
#define TWO_POW_MIN_ONE_THIRD 0.7937005259840998

std::vector<double> Force::repulsion(
	std::vector<double> r_ij,
	double rSquared,
	double radiiSquared,
	double strength,
	std::string typeI,
	std::string typeJ)
{
	if ( (typeI == "soft") && (typeJ == "soft") ) {
		return this->softcoreForce(r_ij, rSquared, radiiSquared, strength);
	}
	if ( (typeI == "lj") && (typeJ == "lj") ) {
		double correctedSigmaSquared = TWO_POW_MIN_ONE_THIRD * radiiSquared;
		return this->LJ1206(r_ij, rSquared, correctedSigmaSquared, strength);
	}
	if ( (typeI == "lj") && (typeJ == "soft") ) {
		double correctedSigmaSquared = TWO_POW_MIN_ONE_THIRD * radiiSquared;
		return this->LJ1206(r_ij, rSquared, correctedSigmaSquared, strength);
	}
	if ( (typeJ == "soft") && (typeJ == "lj") ) {
		double correctedSigmaSquared = TWO_POW_MIN_ONE_THIRD * radiiSquared;
		return this->LJ1206(r_ij, rSquared, correctedSigmaSquared, strength);
	}
	else {
		throw "Particle types are neither 'soft' nor 'lj' !";
	}
}

std::vector<double> Force::softcoreForce(
	std::vector<double> r_ij,
	double rSquared, 
	double radiiSquared,
	double strength)
{
	if ( rSquared > radiiSquared ) {
		std::vector<double> zero = {0.,0.,0.};
		return zero;
	}
	double preFactor = strength * (1. - sqrt(radiiSquared/rSquared));
	std::vector<double> force = {0.,0.,0.};
	force.shrink_to_fit();
	force[0] = preFactor * r_ij[0];
	force[1] = preFactor * r_ij[1]; 
	force[2] = preFactor * r_ij[2];
	return force;
}

std::vector<double> Force::LJ1206(
	std::vector<double> r_ij,
	double rSquared,
	double sigmaSquared,
	double strength)
{
	if ( rSquared > (6.25*sigmaSquared) ) {
		std::vector<double> zero = {0.,0.,0.};
		return zero;
	}
	double preFactor = -4. * strength;
	preFactor *= pow(sigmaSquared/rSquared,6.) * (12./rSquared) 
	           - pow(sigmaSquared/rSquared,3.) * (6. /rSquared);
	std::vector<double> force = {0.,0.,0.};
	force.shrink_to_fit();
	force[0] = preFactor * r_ij[0];
	force[1] = preFactor * r_ij[1];
	force[2] = preFactor * r_ij[2];
	return force;
}
