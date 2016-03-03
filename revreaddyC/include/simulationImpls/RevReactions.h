/* RevReactions.h - child of SimulationImpl that has the hardcoded procedure
 * of performing standard Brownian Dynamics but reactions with an acceptance
 * step. Thus only the run() method is altered compared to the default 
 * implementation. */
#ifndef __REVREACTIONS_H_INCLUDED__
#define __REVREACTIONS_H_INCLUDED__
#include "SimulationImpl.h"

class RevReactions : public SimulationImpl {
public:
 	/* Only the RevReactions constructor is needed.
 	 * For destruction only the SimulationImpl destructor is called
 	 * because it would do the same as the individual destructors
 	 * anyway. */
	RevReactions(World * inWorld, Config * inConfig);
	void run(const unsigned long maxTime);
};
#endif //__REVREACTIONS_H_INCLUDED__