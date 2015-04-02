/* Geometry.cpp */

#include "Geometry.h"

void Geometry::forceEnergy(
	std::vector<double>& force,
	double& energy,
	std::vector<double>& particlePosition,
	double& particleRadius)
{

}

bool Geometry::doesInteract(unsigned int& particleTypeId)
{
	unsigned int max = this->particleTypeIds.size() - 1;
	unsigned int min = 0;
	unsigned int mid = 0;
	while (max >= min) {
		mid = min + (max - min) / 2;
		if (this->particleTypeIds[mid] == particleTypeId) {return true;}
		else if (this->particleTypeIds[mid] < particleTypeId) {min = mid + 1;}
		else {max = mid - 1;}
	}
	return false;
}
