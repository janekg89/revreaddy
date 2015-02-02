/* Force.h
 * author: Christoph Froehner
 * 
 * Handler for calculation of forces and energies.
 */

#ifndef __FORCE_H_INCLUDED__
#define __FORCE_H_INCLUDED__
#include <vector>
#include <math.h>
#include <cmath>
#include <iostream>

class Force
{
	public:
		void repulsionForce(
			std::vector<double>& forceI,//out
			std::vector<double>& r_ij,//in
			double& rSquared,
			double& radiiSquared,
			double& strength,
			unsigned int& typeI,
			unsigned int& typeJ);
		void softcoreForce(
			std::vector<double>& forceI,//out
			std::vector<double>& r_ij,//in
			double& rSquared, 
			double& radiiSquared,
			double& strength);
		void LJ1206Force(
			std::vector<double>& forceI,//out
			std::vector<double>& r_ij,//in
			double& rSquared,
			double& sigmaSquared,
			double& strength);

		void repulsionEnergy(
			double& energy,//out 
			double& rSquared,//in
			double& radiiSquared,
			double& strength,
			unsigned int& typeI,
			unsigned int& typeJ);
		void softcoreEnergy(
			double& energy,//out
			double& rSquared,//in
			double& radiiSquared,
			double& strength);
		void LJ1206Energy(
			double& energy,//out
			double& rSquared,//in
			double& sigmaSquared,
			double& strength);

		void repulsionForceEnergy(
			std::vector<double>& forceI,//out
			double& energy,
			std::vector<double>& r_ij,//in
			double& rSquared,
			double& radiiSquared,
			double& strength,
			unsigned int& typeI,
			unsigned int& typeJ);
		void softcoreForceEnergy(
			std::vector<double>& forceI,//out
			double& energy,
			std::vector<double>& r_ij,//in
			double& rSquared, 
			double& radiiSquared,
			double& strength);
		void LJ1206ForceEnergy(
			std::vector<double>& forceI,//out
			double& energy,
			std::vector<double>& r_ij,//in
			double& rSquared,
			double& sigmaSquared,
			double& strength);
};
#endif // __FORCE_H_INCLUDED__
