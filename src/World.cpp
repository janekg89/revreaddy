/* World.cpp */

#include "World.h"

World::World()
{
	this->cumulativeRuntime     = 0.;
	this->energy                = 0.;
	this->oldEnergy             = 0.;	
	this->acceptProbDynamics    = 1.;
	this->acceptProbReactions   = 1.;
	this->acceptionsDynamics    = 0;
	this->rejectionsDynamics    = 0;
	this->uniqueIdCounter       = 0;
}