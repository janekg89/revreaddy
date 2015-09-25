/* UnimolecularCandidate.h
 * This structure tells that one particleType can undergo
 * one certain reaction, given its reactionId and
 * direction (forwardOrBackward) */

 #ifndef __UNIMOLECULARCANDIDATE_H_INCLUDED__
 #define __UNIMOLECULARCANDIDATE_H_INCLUDED__

 struct UnimolecularCandidate
 {
	public:
	UnimolecularCandidate(
		unsigned inParticleTypeId,
		unsigned inReactionId,
		bool inForwardOrBackward);

	unsigned particleTypeId;
	unsigned reactionId;
	bool forwardOrBackward;
 };

 #endif // __UNIMOLECULARCANDIDATE_H_INCLUDED__