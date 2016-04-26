/* Acceptance.h
 * Record the acceptance probability as a function of time. */

#ifndef __ACCEPTANCE_H_INCLUDED__
#define __ACCEPTANCE_H_INCLUDED__
#include <vector>
#include <string>
#include "World.h"
#include "Config.h"
#include "Observable.h"

// TODO make two Acceptances one for Reactions, one for Diffusion
class Acceptance : public Observable
{
public:
	Acceptance(
		unsigned int inRecPeriod,
		unsigned int inClearPeriod,
		std::string inFilename,
		bool inReactionsOrDiffusion);
	void configure(World * world, Config * config);
	void record(World * world, double t);
	void writeToH5();
	void writeToDat();
private:
	std::vector<double> acceptanceProbs;
	std::vector<double> times;
	bool reactionsOrDiffusion;
};

#endif // __ACCEPTANCE_H_INCLUDED__