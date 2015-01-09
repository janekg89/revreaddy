/* utils.cpp */

#include "utils.h"

void printVector(std::vector<double> vec)
{
	std::cout << "[ ";
	for (int i = 0; i < vec.size(); i++)
	{
		std::cout << vec[i] << " ";
	}
	std::cout << "]" << std::endl;
}

double squaredDistance(std::array<double, 3> arr1, std::array<double, 3> arr2)
{
	double result;
	result = ( arr1[0] - arr2[0] ) * ( arr1[0] - arr2[0] );
	result += ( arr1[1] - arr2[1] ) * ( arr1[1] - arr2[1] );
	result += ( arr1[2] - arr2[2] ) * ( arr1[2] - arr2[2] );
	return result;
}

// Return the minimum distance vector, pointing from r_i to r_j
std::array<double,3>
getMinDistanceVector(
	std::array<double,3> r_i,
	std::array<double,3> r_j,
	bool isPeriodic,
	double boxsize)
{
	double dx, dy, dz;
	std::array<double,3> r_ij;
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
	return r_ij;
}
