/* TypeDict.cpp */

#include "TypeDict.h"

TypeDict::TypeDict()
{
	this->numberOfTypes = 0;
}

TypeDict::~TypeDict()
{

}

void TypeDict::newType(string name, double radius, double diffusionConst)
{
	this->names.push_back(name);
	this->radii.push_back(radius);
	this->diffusionConstants.push_back(diffusionConst);
	this->numberOfTypes += 1;
}

unsigned int TypeDict::getNumberOfTypes()
{
	return this->numberOfTypes;
}
