/* Geometry.h */

#ifndef __GEOMETRY_H_INCLUDED__
#define __GEOMETRY_H_INCLUDED__
#include <vector>
#include "Particle.h"

class Geometry
{
	public:
		/* SORTED list of particle type ids that interact
		 * with geometry. */
		std::vector<unsigned int> particleTypeIds;

		void forceEnergy(
			std::vector<double>& force,
			double& energy,
			std::vector<double>& particlePosition,
			double& particleRadius);
		/* binary search within particleTypeIds */
		bool doesInteract(unsigned int& particleTypeId);
};

#endif
