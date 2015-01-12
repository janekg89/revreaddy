/* Force.h
 * author: Christoph Froehner
 * 
 * Handler for calculation of forces.
 */

#ifndef __FORCE_H_INCLUDED__
#define __FORCE_H_INCLUDED__
#include <array>
#include <math.h>
#include <cmath>
#include <string>

class Force
{
	public:
		std::array<double,3> repulsion(
			std::array<double,3> r_ij,
			double rSquared,
			double radiiSquared,
			double strength,
			std::string typeI,
			std::string typeJ);
		std::array<double,3> softcoreForce(
			std::array<double,3> r_ij,
			double rSquared, 
			double radiiSquared,
			double strength);
		std::array<double,3> LJ1206(
			std::array<double,3> r_ij,
			double rSquared,
			double sigmaSquared,
			double strength);
};
#endif // __FORCE_H_INCLUDED__
