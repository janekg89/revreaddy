/* Potential.cpp
 * author: Christoph Froehner
 */

#include "Potential.h"

// TODO implement force correctly
std::array<double, 3> Potential::softcore(std::array<double, 3> r_ij, double rSquared, 
double cutoffSquared, double strength)
{
	std::array<double, 3> force = {0.,0.,0.};
	return force;
}
