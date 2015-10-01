/* ParticleInteraction.h
 * 
 * This class will have children that implement a special
 * type of force, e.g. soft-repulsion, soft-attraction
 * for ONE pair of particle types. */

#ifndef __PARTICLEINTERACTION_H_INCLUDED__
#define __PARTICLEINTERACTION_H_INCLUDED__
#include <string>
#include <vector>
#include <iostream>

class ParticleInteraction
{
	public:
		ParticleInteraction();
		virtual ~ParticleInteraction();

		std::string name;
		std::string type; //"softRepulsion", "softAttraction", "LennardJones"
		/* affectedTuple is always of length 2 and always sorted */
		std::vector<unsigned int> affectedTuple; 
		std::vector<double> parameters;
		/* cutoff is important only for neighborlattice construction */
		double cutoff;
		/* reference counting for usage with reactions to 
		 * keep track if the ParticleInteraction is still 
		 * referenced by some Reaction */
		unsigned references;
		
		bool isAffected(unsigned int i, unsigned int j);
		virtual void calculateForceEnergy(
			std::vector<double>& forceI,//out
			double& energy,
			std::vector<double>& r_ij,//in
			double& rSquared,
			double& radiiSquared);
		virtual double calculateEnergy(
			double rSquared, // in
			double radiiSquared); //in
		void incrementRef();
		void decrementRef();
};

#endif // __PARTICLEINTERACTION_H_INCLUDED__