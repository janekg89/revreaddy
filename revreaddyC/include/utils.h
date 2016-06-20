/* utils.h */

#ifndef __UTILS_H_INCLUDED__
#define __UTILS_H_INCLUDED__
#include <vector>
#include <iostream>
#include <cmath>

class Utils {
public:
	Utils()
	{
		dx = 0;
		dy = 0;
		dz = 0;
	}
	double dx, dy, dz;
	void getMinDistanceVector(
		std::vector<double>& r_ij,
		std::vector<double>& r_i,
		std::vector<double>& r_j,
		bool& isPeriodic,
		double& boxsize);
	void getMinDistanceSquared(
		double& distance,
		std::vector<double>& r_i,
		std::vector<double>& r_j,
		bool& isPeriodic,
		double& boxsize);
};

inline void Utils::getMinDistanceVector(
	std::vector<double>& r_ij,
	std::vector<double>& r_i,
	std::vector<double>& r_j,
	bool& isPeriodic,
	double& boxsize
)
{
	dx = r_j[0] - r_i[0];
	dy = r_j[1] - r_i[1];
	dz = r_j[2] - r_i[2];
	if ( isPeriodic )
	{
		// copysign returns a value with the magnitude of the first arg and
		// the sign of the second arg
		if ( fabs(dx) > 0.5*boxsize ) { dx -= copysign(boxsize, dx); }
		if ( fabs(dy) > 0.5*boxsize ) { dy -= copysign(boxsize, dy); }
		if ( fabs(dz) > 0.5*boxsize ) { dz -= copysign(boxsize, dz); }
	}
	r_ij[0] = dx;
	r_ij[1] = dy;
	r_ij[2] = dz;
}

inline void Utils::getMinDistanceSquared(
	double& distance,
	std::vector<double>& r_i,
	std::vector<double>& r_j,
	bool& isPeriodic,
	double& boxsize)
{
	dx = r_j[0] - r_i[0];
	dy = r_j[1] - r_i[1];
	dz = r_j[2] - r_i[2];
	if ( isPeriodic )
	{
		// copysign returns a value with the magnitude of the first arg and
		// the sign of the second arg
		if ( fabs(dx) > 0.5*boxsize ) { dx -= copysign(boxsize, dx); }
		if ( fabs(dy) > 0.5*boxsize ) { dy -= copysign(boxsize, dy); }
		if ( fabs(dz) > 0.5*boxsize ) { dz -= copysign(boxsize, dz); }
	}
	distance = dx*dx + dy*dy + dz*dz;
}

#endif //__UTILS_H_INCLUDED__