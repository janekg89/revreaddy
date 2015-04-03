/* DoubleWellZ.cpp */

#include "DoubleWellZ.h"
#define SHIFT 0.1125 // the term 3a/8 if a = 0.3
#define SCALEFACTOR 1.472667648

DoubleWellZ::DoubleWellZ(double InDistanceMinima, double InStrength)
{
	this->distanceMinima = InDistanceMinima;
	this->strength = InStrength;
	this->scale = this->distanceMinima / SCALEFACTOR;
}

void DoubleWellZ::forceEnergy(
	std::vector<double>& force,
	double& energy,
	std::vector<double>& particlePosition,
	double& particleRadius)
{
	double term = particlePosition[2] / this->scale - SHIFT;
	double energyTerm = 0.;
	double forceTerm = -2. * term;
	term       *= term;
	forceTerm  += 0.9 * term;
	energyTerm  = -1. * term;
	term       *= term;
	forceTerm  += 4.  * term;
	energyTerm += 0.3 * term;
	term       *= term;
	energyTerm += term;
	forceTerm  *= -1. * strength;
	energy = energyTerm * strength;
	force[0] = 0.;
	force[1] = 0.;
	force[2] = forceTerm;
}
