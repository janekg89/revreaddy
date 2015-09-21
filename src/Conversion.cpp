/* Conversion.cpp */

#include "Conversion.h"

Conversion::Conversion(
	std::string inName,
	std::vector<unsigned int> inForwardTypes,
	std::vector<unsigned int> inBackwardTypes,
	double inForwardRate,
	double inBackwardRate)
{
	this->name = inName;
	this->forwardTypes = inForwardTypes;
	this->backwardTypes = inBackwardTypes;
	this->forwardRate = inForwardRate;
	this->backwardRate = inBackwardRate;
}

Conversion::~Conversion()
{

}

double Conversion::performForward(
	std::vector<unsigned long int> particleIndices,
	World * world)
{
	return 1.;
}

double Conversion::performBackward(
	std::vector<unsigned long int> particleIndices,
	World * world)
{
	return 1.;
}