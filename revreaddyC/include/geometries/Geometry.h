/* Geometry.h */

#ifndef __GEOMETRY_H_INCLUDED__
#define __GEOMETRY_H_INCLUDED__
#include <algorithm>
#include <vector>
#include "Particle.h"

class Geometry
{
	public:
		Geometry();
		virtual ~Geometry();
		/* SORTED list of particle type ids that interact
		 * with geometry. */
		std::vector<unsigned int> particleTypeIds;

		/* the actual implementation depends on the children of Geometry */
		virtual void forceEnergy(
			std::vector<double>& force,
			double& energy,
			std::vector<double>& particlePosition,
			double& particleRadius);
		/* binary search within particleTypeIds */
		bool doesInteract(unsigned int& particleTypeId);
};

#endif
