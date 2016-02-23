/* ReactionEvent.cpp */

#include "ReactionEvent.h"

ReactionEvent::ReactionEvent(
	unsigned int inReactionId,
	bool inForwardOrBackward,
	std::vector<unsigned long long> inParticipants)
{
	this->reactionId = inReactionId;
	this->forwardOrBackward = inForwardOrBackward;
	this->participants = inParticipants;
}
