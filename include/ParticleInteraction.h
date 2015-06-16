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
		ParticleInteraction();
		virtual ~ParticleInteraction();

		std::string name;
		std::string type; //"softRepulsion", "softAttraction", "LennardJones"
		/* affectedTuple is always of length 2 and always sorted */
		std::vector<unsigned int> affectedTuple; 
		double cutoff;
		
		bool isAffected(unsigned int i, unsigned int j);
		virtual void calculateForceEnergy(
			std::vector<double>& forceI,//out
			double& energy,
			std::vector<double>& r_ij,//in
			double& rSquared,
			double& radiiSquared);
};

#endif // __PARTICLEINTERACTION_H_INCLUDED__
