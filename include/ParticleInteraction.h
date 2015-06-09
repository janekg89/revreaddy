/* ParticleInteraction.h
 * 
 * This class will have children that implement a special
 * type of force, e.g. soft-repulsion, soft-attraction
 * for ONE pair of particle types. */

#ifndef __PARTICLEINTERACTION_H_INCLUDED__
#define __PARTICLEINTERACTION_H_INCLUDED__
#include <string>
#include <vector>

class ParticleInteraction
{
	public:
		std::string name;
		std::string type;
		std::vector<unsigned int> affectedTuple;
		
		bool isAffected(unsigned int i, unsigned int j);
		void calculateForceEnergy(
			std::vector<double>& forceI,//out
			double& energy,
			std::vector<double>& r_ij,//in
			double& rSquared,
			double& radiiSquared,
			double& strength);
}

#endif // __PARTICLEINTERACTION_H_INCLUDED__
