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

inline void print3DArray(std::array<double, 3> arr)
{
	std::cout << "[ ";
	for (int i = 0; i < arr.size(); i++)
	{
		std::cout << arr[i] << " ";
	}
	std::cout << "]" << std::endl;
}

inline void printVector(std::vector<double> vec)
{
	std::cout << "[ ";
	for (int i = 0; i < vec.size(); i++)
	{
		std::cout << vec[i] << " ";
	}
	std::cout << "]" << std::endl;
}

inline double squaredDistance(std::array<double, 3> arr1, std::array<double, 3> arr2)
{
	double result;
	result = ( arr1[0] - arr2[0] ) * ( arr1[0] - arr2[0] );
	result += ( arr1[1] - arr2[1] ) * ( arr1[1] - arr2[1] );
	result += ( arr1[2] - arr2[2] ) * ( arr1[2] - arr2[2] );
	return result;
}

#endif
