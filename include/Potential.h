/* Potential.h
 * author: Christoph Froehner
 * 
 * Handler for calculation of forces and energies.
 */

#ifndef __POTENTIAL_H_INCLUDED__
#define __POTENTIAL_H_INCLUDED__
#include <array>
#include <math.h>
#include <cmath>
#include <string>

class Potential
{
	public:
		std::array<double,3> repulsion(std::array<double,3> r_ij, double rSquared,
			double radiiSquared, double strength, std::string typeI, std::string typeJ);
		std::array<double,3> softcoreForce(std::array<double, 3> r_ij, double rSquared, 
			double radiiSquared, double strength);
		std::array<double,3> LJ1206(std::array<double,3> r_ij, double rSquared,
			double sigmaSquared, double strength);
};
#endif
