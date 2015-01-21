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
	std::vector<double> arr1,
	std::vector<double> arr2
);

void getMinDistanceVector(
	std::vector<double>& r_ij,
	std::vector<double>& r_i,
	std::vector<double>& r_j,
	bool& isPeriodic,
	double& boxsize
);

#endif
