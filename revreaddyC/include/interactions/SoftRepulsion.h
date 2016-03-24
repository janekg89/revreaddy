/* SoftRepulsion.h 
 * 
 * Implement the repulsive harmonic force as a child
 * of the ParticleInteraction class. */

#ifndef __SOFTREPULSION_H_INCLUDED__
#define __SOFTREPULSION_H_INCLUDED__
#include <math.h>
#include <cmath>
#include <algorithm>
#include "ParticleInteraction.h"
#include "logging.h"

class SoftRepulsion : public ParticleInteraction {
public:
	SoftRepulsion(std::string inName, std::array<unsigned,2> inAffectedTuple, double inRepulsionStrength);
	
	void calculateForceEnergy(
		std::vector<double>& forceI, //out
		double& energy,
		const std::vector<double>& r_ij, //in
		const double& rSquared,
		const double& radiiSquared);
	double calculateEnergy(
		const double rSquared, // in
		const double radiiSquared); //in
	/* members that are once allocated to avoid
	 * on-the-run allocation */
	double c;
	double a;
};

#endif // __SOFTREPULSION_H_INCLUDED__