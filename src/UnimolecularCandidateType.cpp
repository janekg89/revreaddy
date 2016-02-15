/* UnimolecularCandidate.cpp */

#include "UnimolecularCandidate.h"

UnimolecularCandidate::UnimolecularCandidate(
	unsigned inParticleTypeId,
	unsigned inReactionId,
	bool inForwardOrBackward)
{
	this->particleTypeId = inParticleTypeId;
	this->reactionId = inReactionId;
	this->forwardOrBackward = inForwardOrBackward;
}