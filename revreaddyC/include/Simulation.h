/* Simulation.h
 * author: Christoph Froehner
 *
 * This is the main class, which uses all other classes.
 * Simulation performs the time loop which propagates
 * the particles and also holds information that concern
 * a particular realization, e.g. the list of observables.
 * The system dependent variables are stored in World and
 * Config. The implementation is found in SimulationImpl. */

#ifndef __SIMULATION_H_INCLUDED__
#define __SIMULATION_H_INCLUDED__
#include <string>
#include <vector>
#include <array>
#include "Config.h"
#include "World.h"
#include "SimulationImpl.h"
#include "RevReactions.h"
#include "RevDiffusion.h"
#include "RevDiffusionRevReactions.h"
#include "FractionalDiffusion.h"
#include "logging.h"

class Simulation {
public:
	Simulation(World * inWorld, Config * inConfig, std::string whichImpl);
	~Simulation();
	/* Start the simulation. Iterate for maxTime timesteps.*/
	void run(const unsigned long maxTime);
	bool getUseNeighborlist();
	void setUseNeighborlist(const bool inUseNeighborlist);

	/* Observable methods */
	void writeAllObservablesToFile();
	std::string showObservables();
	void deleteAllObservables();
	void new_Trajectory(unsigned long recPeriod, std::string filename);
	void new_TrajectoryUnique(unsigned long recPeriod, unsigned long clearPeriod, std::string filename);
	void new_RadialDistribution(unsigned long recPeriod, std::string filename, std::vector<double> ranges, std::vector< std::array<unsigned,2> > considered, std::vector<unsigned long> recordingRange);
	void new_MeanSquaredDisplacement(unsigned long recPeriod, std::string filename,	unsigned particleTypeId);
	void new_ProbabilityDensity(unsigned long recPeriod, std::string filename, unsigned pTypeId, std::vector<double> range, unsigned int coord);
	void new_Energy(unsigned long recPeriod, std::string filename);
	void new_Acceptance(unsigned long recPeriod, std::string filename, bool reactionsOrDiffusion);
	void new_ParticleNumbers(unsigned long recPeriod, std::string filename,	unsigned particleTypeId);
	void new_Increments(unsigned long recPeriod, unsigned long clearPeriod, std::string filename, unsigned particleTypeId);
	void new_ReactionCounter(unsigned long recPeriod, std::string filename, std::string reactionName);

private:
	SimulationImpl * impl;
};

#endif // __SIMULATION_H_INCLUDED__