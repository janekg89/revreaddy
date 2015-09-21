/* Energy.h
 * 
 * Record the total energy of the system as a function of time.
 * This observable has a special role, because it is granted
 * knowledge of the Simulation object. This is necessary to
 * obtain the energy variable from Simulation. */

#ifndef __ENERGY_H_INCLUDED__
#define __ENERGY_H_INCLUDED__
#include <vector>
#include <string>
#include "Observable.h"

class Energy : public Observable
{
	public:
		std::vector<double> energies;
		std::vector<double> times;

		Energy(
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


#endif // __ENERGY_H_INCLUDED__
