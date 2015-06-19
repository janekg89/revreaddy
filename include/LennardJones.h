/* LennardJones.h
 *
 * LennardJones interaction as a child of the
 * ParticleInteraction class. */

#ifndef __LENNARDJONES_H_INCLUDED__
#define __LENNARDJONES_H_INCLUDED__
#include <math.h>
#include <cmath>
#include "ParticleInteraction.h"

class LennardJones : public ParticleInteraction
{
	public:
		LennardJones(
			std::string inName,
			std::vector<unsigned int> inAffectedTuple,
			double inEpsilon);

		void calculateForceEnergy(
			std::vector<double>& forceI, //out
			double& energy,
			std::vector<double>& r_ij, //in
			double& rSquared,
			double& radiiSquared);
};
		
#endif // __LENNARDJONES_H_INCLUDED__
