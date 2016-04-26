/* UnimolecularCandidateType.h
 * This structure tells that one particleType can undergo
 * one certain reaction, given its reactionId and
 * direction (forwardOrBackward) */

 #ifndef __UNIMOLECULARCANDIDATETYPE_H_INCLUDED__
 #define __UNIMOLECULARCANDIDATETYPE_H_INCLUDED__

 struct UnimolecularCandidateType {
	UnimolecularCandidateType(unsigned pTypeId,	unsigned reactId, bool forwBackw) 
		: particleTypeId(pTypeId), reactionId(reactId), forwardOrBackward(forwBackw) { }
	unsigned particleTypeId;
	unsigned reactionId;
	bool forwardOrBackward;
 };

 #endif // __UNIMOLECULARCANDIDATETYPE_H_INCLUDED__