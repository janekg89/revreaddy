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
	std::array<double, 3> arr1,
	std::array<double, 3> arr2
);
std::array<double,3> getMinDistanceVector(
	std::array<double,3> r_i,
	std::array<double,3> r_j,
	bool isPeriodic,
	double boxsize
);

#endif
