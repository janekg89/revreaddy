/* Potential.h
 * author: Christoph Froehner
 * 
 * Handler for calculation of forces and energies.
 */

#ifndef __POTENTIAL_H_INCLUDED__
#define __POTENTIAL_H_INCLUDED__
#include <array>

// TODO How to handle calculation of forces and energies, while calculating the distance
//  only once.
class Potential
{
	public:
		std::array<double,3> softcore(std::array<double, 3> r_ij, double rSquared, 
		double cutoffSquared, double strength);
};
#endif
