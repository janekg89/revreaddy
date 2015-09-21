/* Acceptance.h
 *
 * Record the acceptance probability as a function of time.
 * This observable has a special role, because it is granted
 * knowledge of the Simulation object. This is necessay to
 * obtain the currentAcceptance variable from Simulation. */

#ifndef __ACCEPTANCE_H_INCLUDED__
#define __ACCEPTANCE_H_INCLUDED__
#include <vector>
#include <string>
#include "Observable.h"

// TODO make two Acceptances one for Reactions, one for Dynamics
class Acceptance : public Observable
{
	public:
		std::vector<double> acceptanceProbs;
		std::vector<double> times;
		
		Acceptance(
			unsigned int inRecPeriod,
			unsigned int inClearPeriod,
			std::string inFilename);

		void record(
			World * world,
			double t);
		void writeBufferToFile();
		void writeBufferToH5();
		void writeBufferToDat();
};

#endif // __ACCEPTANCE_H_INCLUDED__
