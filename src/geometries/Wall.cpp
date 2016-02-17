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
	this->R = 0.;
	this->preFactor = 0.;
}

void Wall::forceEnergy(
	std::vector<double>& force,
	double& energy,
	std::vector<double>& particlePosition,
	double& particleRadius)
{
	this->R = ( particlePosition[0] - point[0] ) * normal[0]
	         + ( particlePosition[1] - point[1] ) * normal[1]
	         + ( particlePosition[2] - point[2] ) * normal[2];
	//         - particleRadius;
	if (this->R > 0.) {
		force[0] = 0.;
		force[1] = 0.;
		force[2] = 0.;
		energy = 0.;
		return;
	}
	else {
		this->preFactor = strength * this->R;
		energy = this->preFactor * this->R;
		this->preFactor *= -2.;
		force[0] = this->preFactor * normal[0];
		force[1] = this->preFactor * normal[1];
		force[2] = this->preFactor * normal[2];
	}
}