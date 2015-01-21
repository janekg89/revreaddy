/* Force.h
 * author: Christoph Froehner
 * 
 * Handler for calculation of forces.
 */

#ifndef __FORCE_H_INCLUDED__
#define __FORCE_H_INCLUDED__
#include <vector>
#include <math.h>
#include <cmath>
#include <string>
#include <iostream>

class Force
{
	public:
		void repulsionForce(
			std::vector<double>& forceI,
			std::vector<double>& r_ij,
			double& rSquared,
			double& radiiSquared,
			double& strength,
			std::string& typeI,
			std::string& typeJ);
		void softcoreForce(
			std::vector<double>& forceI,
			std::vector<double>& r_ij,
			double& rSquared, 
			double& radiiSquared,
			double& strength);
		void LJ1206Force(
			std::vector<double>& forceI,
			std::vector<double>& r_ij,
			double& rSquared,
			double& sigmaSquared,
			double& strength);
		double repulsionEnergy(
			double rSquared,
			double radiiSquared,
			double strength,
			std::string typeI,
			std::string typeJ);
		double softcoreEnergy(
			double rSquared,
			double radiiSquared,
			double strength);
		double LJ1206Energy(
			double rSquared,
			double sigmaSquared,
			double strength);
};
#endif // __FORCE_H_INCLUDED__
