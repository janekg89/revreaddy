/* utils.cpp */

#include "utils.h"

void printVector(std::vector<double> vec)
{
	std::cout << "[ ";
	for (unsigned i = 0; i < vec.size(); i++) {
		std::cout << vec[i] << " ";
	}
	std::cout << "]" << std::endl;
}