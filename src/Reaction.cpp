/* Reaction.cpp */

#include "Reaction.h"

Reaction::Reaction()
{

}

Reaction::~Reaction()
{
	
}

// TODO could be faster, when types are checked for inequality first ???
// TODO also it looks ugly!
bool Reaction::isAffectedForward(std::vector<unsigned int> types)
{
	if ( this->forwardTypes.size() != types.size() ) {return false;}
	if ( types.size()==2 ) {
		if ( (types[0]==this->forwardTypes[0]) 
		     && (types[1]==this->forwardTypes[1]) ) {
			return true;
		}
		else if ( (types[0]==this->forwardTypes[1]) 
		          && (types[1]==this->forwardTypes[0]) ) {
			return true;
		}
		else {return false;}
	}
	else if ( types.size()==1 ) {return (types[0]==this->forwardTypes[0]);}
	else {return false;}
}

bool Reaction::isAffectedForward(unsigned int type)
{
	if ( this->forwardTypes.size()!=1 ) {return false;}
	else {return (type==this->forwardTypes[0]);}
}

bool Reaction::isAffectedBackward(std::vector<unsigned int> types)
{
	if ( this->backwardTypes.size() != types.size() ) {return false;}
	if ( types.size()==2 ) {
		if ( (types[0]==this->backwardTypes[0]) 
		     && (types[1]==this->backwardTypes[1]) ) {
			return true;
		}
		else if ( (types[0]==this->backwardTypes[1]) 
		          && (types[1]==this->backwardTypes[0]) ) {
			return true;
		}
		else {return false;}
	}
	else if ( types.size()==1 ) {return (types[0]==this->backwardTypes[0]);}
	else {return false;}
}

bool Reaction::isAffectedBackward(unsigned int type)
{
	if ( this->backwardTypes.size()!=1 ) {return false;}
	else {return (type==this->backwardTypes[0]);}
}

double Reaction::performForward(
	std::vector<unsigned long int> particleIndices,
	World * world)
{
	return 1.;
}

double Reaction::performBackward(
	std::vector<unsigned long int> particleIndices,
	World * world)
{
	return 1.;
}