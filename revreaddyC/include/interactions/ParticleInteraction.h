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
#include "Exception.h"

class ParticleInteraction {
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
	
	bool isAffected(unsigned int i, unsigned int j);
	virtual void calculateForceEnergy(
		std::vector<double>& forceI,//out
		double& energy,
		const std::vector<double>& r_ij,//in
		const double& rSquared,
		const double& radiiSquared);
	virtual double calculateEnergy(
		const double rSquared, // in
		const double radiiSquared); //in
};

#endif // __PARTICLEINTERACTION_H_INCLUDED__