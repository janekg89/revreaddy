/* TypeDict.cpp */

#include "TypeDict.h"

TypeDict::TypeDict()
{
	this->numberOfTypes = 0;
}

TypeDict::~TypeDict()
{

}

void TypeDict::newType(
	std::string name,
	double radius,
	double diffusionConst,
	double reactionRadius,
	unsigned int forceType)
{
	this->names.push_back(name);
	this->radii.push_back(radius);
	this->diffusionConstants.push_back(diffusionConst);
	this->reactionRadii.push_back(reactionRadius);
	this->forceTypes.push_back(forceType);
	this->numberOfTypes += 1;
}

unsigned int TypeDict::getNumberOfTypes()
{
	return this->numberOfTypes;
}
