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
#include "Simulation.h"

class Acceptance : public Observable
{
	public:
		std::vector<double> acceptanceProbs;
		std::vector<double> times;
		/* pointer to the simulation, owning this observable */
		Simulation * simulation;
		
		Acceptance(
			unsigned int inRecPeriod,
			unsigned int inClearPeriod,
			Simulation * inSimulation,
			std::string inFilename);

		void record(
			std::vector<Particle>& activeParticles,
			double t);
		void writeBufferToFile();
		void writeBufferToH5();
		void writeBufferToDat();
};

#endif // __ACCEPTANCE_H_INCLUDED__
