/* Geometry.cpp */
#include "Geometry.h"
Geometry::Geometry() {}
Geometry::~Geometry() {}

void Geometry::forceEnergy(
	std::vector<double>& force,
	double& energy,
	std::vector<double>& particlePosition,
	double& particleRadius)
{}

bool Geometry::doesInteract(unsigned int& particleTypeId) {
	// as far as particleTypeIds is not sorted, we do a linear search here
	auto result = std::find(this->particleTypeIds.begin(), this->particleTypeIds.end(), particleTypeId);
	return result != particleTypeIds.end();
	// max min mid have to be signed, otherwise max might underflow undetected
	/*signed int max = this->particleTypeIds.size() - 1;
	signed int min = 0;
	signed int mid = 0;
	while (max >= min) {
		mid = min + (max - min) / 2;
		if (this->particleTypeIds[mid] == particleTypeId) {return true;}
		else if (this->particleTypeIds[mid] < particleTypeId) {min = mid + 1;}
		else {max = mid - 1;}
	}
	return false;*/
}