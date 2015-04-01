/* Geometry.h */

#ifndef __GEOMETRY_H_INCLUDED__
#define __GEOMETRY_H_INCLUDED__
#include <vector>
#include "Particle.h"

class Geometry
{
	public:
		forceEnergy(
			std::vector<double>& force,
			double& energy,
			std::vector<double>& particlePosition,
			double& particleRadius);
};

#endif
