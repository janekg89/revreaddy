/* utils.h
 * author: Christoph Froehner
 *
 * several functions needed everywhere.
 */

#ifndef __UTILS_H_INCLUDED__
#define __UTILS_H_INCLUDED__
#include <array>
#include <vector>
#include <iostream>
#include <cmath>

template<std::size_t SIZE>
void printArray(std::array<double, SIZE>& arr)
{
	std::cout << "[ ";
	for(auto& entry : arr)
	{
		std::cout << entry << " ";
	}
	std::cout << "]" << std::endl;
}

void printVector(std::vector<double> vec);

double squaredDistance(
	std::vector<double>& arr1,
	std::vector<double>& arr2
);

inline void getMinDistanceVector(
	std::vector<double>& r_ij,
	std::vector<double>& r_i,
	std::vector<double>& r_j,
	bool& isPeriodic,
	double& boxsize
)
{
    double dx, dy, dz;
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

#endif
