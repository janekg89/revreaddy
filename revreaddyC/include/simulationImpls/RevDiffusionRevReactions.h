/* RevDiffusionRevReaction.h - child of SimulationImpl that performs both
 * the diffusion and reactions reversibly. I.e. with an acceptance step. */
#ifndef __REVDIFFUSIONREVREACTIONS_H_INCLUDED__
#define __REVDIFFUSIONREVREACTIONS_H_INCLUDED__
#include "SimulationImpl.h"

class RevDiffusionRevReactions : public SimulationImpl {
public:
	RevDiffusionRevReactions(World * inWorld, Config * inConfig);
	void run(const unsigned long maxTime);
};
#endif//__REVDIFFUSIONREVREACTIONS_H_INCLUDED__