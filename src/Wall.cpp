/* Wall.cpp */

#include "Wall.h"

Wall::Wall(
	std::vector<double>& InNormal,
	std::vector<double>& InPoint,
	double& InStrength,
	std::vector<unsigned int>& InParticleTypeIds)
{
	this->normal = InNormal;
	this->point = InPoint;
	this->strength = InStrength;
	this->particleTypeIds = InParticleTypeIds;
}

void Wall::forceEnergy(
	std::vector<double>& force,
	double& energy,
	std::vector<double>& particlePosition,
	double& particleRadius)
{
	double R = ( particlePosition[0] - point[0] ) * normal[0]
	         + ( particlePosition[1] - point[1] ) * normal[1]
	         + ( particlePosition[2] - point[2] ) * normal[2]
	         - particleRadius;
	if (R > 0.) {return;}
	else {
		double preFactor = strength * R;
		energy = preFactor * R;
		preFactor *= -2.;
		force[0] = preFactor * normal[0];
		force[1] = preFactor * normal[1];
		force[2] = preFactor * normal[2];
	}
}
