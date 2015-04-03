/* DoubleWellZ.h
 * 
 * Realise a double well potential, depending on z coordinate of the
 * following form:
 *
 * V(z) =  e*( ( z/L - 3a/8 )**4 - ( z/L - 3a/8 )**2 + a * ( z/L - 3a/8 )**3 ) 
 * 
 * where e = strength, L = scale = distanceMinima / 1.47267 if a = 0.3.
 */

#ifndef __DOUBLEWELLZ_H_INCLUDED__
#define __DOUBLEWELLZ_H_INCLUDED__
#include <vector>
#include "Geometry.h"

class DoubleWellZ : public Geometry
{
	public:
		double distanceMinima;
		double scale;
		double strength;

		DoubleWellZ(double InDistanceMinima, double InStrength);

		void forceEnergy(
			std::vector<double>& force, //out
			double& energy, //out
			std::vector<double>& particlePosition, //in
			double& particleRadius); //in
};

#endif // __DOUBLEWELLZ_H_INCLUDED__
