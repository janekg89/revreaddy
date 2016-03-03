/* RevReactions.h - child of SimulationImpl that has the hardcoded procedure
 * of performing standard Brownian Dynamics but reactions with an acceptance
 * step. Thus only the run() method is altered compared to the default 
 * implementation. */
#ifndef __REVREACTIONS_H_INCLUDED__
#define __REVREACTIONS_H_INCLUDED__
#include "SimulationImpl.h"

class RevReactions : public SimulationImpl {
public:
	RevReactions(World * inWorld, Config * inConfig);
	~RevReactions();
	void run(const unsigned long maxTime);
};
#endif //__REVREACTIONS_H_INCLUDED__