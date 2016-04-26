/* Wall.h */

#ifndef __WALL_H_INCLUDED__
#define __WALL_H_INCLUDED__
#include <vector>
#include <iostream>
#include "Geometry.h"

class Wall : public Geometry {
public:
	std::vector<double> normal;
	std::vector<double> point;
	double strength;
	
	Wall(std::string inName, std::vector<double>& InNormal, std::vector<double>& InPoint, double& InStrength, std::vector<unsigned int>& InParticleTypeIds);

	void forceEnergy(
		std::vector<double>& force, //out
		double& energy, //out
		std::vector<double>& particlePosition, //in
		double& particleRadius); // in

private:
	double R; // distance from the wall. calculated during runtime
	double preFactor;
};

#endif //__WALL_H_INCLUDED__