/* LennardJones.h - LennardJones interaction as a child of the ParticleInteraction class. */

#ifndef __LENNARDJONES_H_INCLUDED__
#define __LENNARDJONES_H_INCLUDED__
#include <math.h>
#include <cmath>
#include "ParticleInteraction.h"
#include "logging.h"

class LennardJones : public ParticleInteraction {
public:
	LennardJones(std::string inName, std::vector<unsigned int> inAffectedTuple,	double inEpsilon);

	void calculateForceEnergy(
		std::vector<double>& forceI, //out
		double& energy,
		const std::vector<double>& r_ij, //in
		const double& rSquared,
		const double& radiiSquared);
	double calculateEnergy(
		const double rSquared, // in
		const double radiiSquared); //in
};
		
#endif // __LENNARDJONES_H_INCLUDED__