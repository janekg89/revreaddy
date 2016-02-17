/* Energy.h
 * Record the total energy of the system as a function of time.
 * This observable has a special role, because it is granted
 * knowledge of the Simulation object. This is necessary to
 * obtain the energy variable from Simulation. */

#ifndef __ENERGY_H_INCLUDED__
#define __ENERGY_H_INCLUDED__
#include <vector>
#include <string>
#include "World.h"
#include "Config.h"
#include "Observable.h"

class Energy : public Observable
{
public:
	Energy(unsigned inRecPeriod, unsigned inClearPeriod, std::string inFilename);
	void configure(World * world, Config * config);
	void record(World * world, double t);
	void writeBufferToFile();
	void writeBufferToH5();
	void writeBufferToDat();
private:
	std::vector<double> energies;
	std::vector<double> times;
};

#endif // __ENERGY_H_INCLUDED__