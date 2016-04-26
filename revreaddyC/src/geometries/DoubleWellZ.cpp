/* DoubleWellZ.cpp */

#include "DoubleWellZ.h"
#define SHIFT 0.1125 // the term 3a/8 if a = 0.3
#define SCALEFACTOR 1.472667648

DoubleWellZ::DoubleWellZ(std::string inName, double inDistanceMinima, double inStrength, std::vector<unsigned int>& inParticleTypeIds) {
	this->distanceMinima = inDistanceMinima;
	this->strength = inStrength;
	this->scale = this->distanceMinima / SCALEFACTOR;
	this->particleTypeIds = inParticleTypeIds;
	this->name = inName;
	this->type = "DoubleWellZ";
}

void DoubleWellZ::forceEnergy(
	std::vector<double>& force,
	double& energy,
	std::vector<double>& particlePosition,
	double& particleRadius)
{
	double term = ( particlePosition[2] - SHIFT ) / this->scale;
	double tempTerm = term;
	double energyTerm = 0.;
	double forceTerm = -2. * tempTerm;
	tempTerm   *= term;
	forceTerm  += 0.9 * tempTerm;
	energyTerm  = -1. * tempTerm;
	tempTerm   *= term;
	forceTerm  += 4.  * tempTerm;
	energyTerm += 0.3 * tempTerm;
	tempTerm   *= term;
	energyTerm += tempTerm;
	forceTerm  *= -1. * strength;
	energy = energyTerm * strength;
	force[0] = 0.;
	force[1] = 0.;
	force[2] = forceTerm;
}
