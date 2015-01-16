/* TrajectorySingle.h 
 * Record the position of a single particle.
 * This observable is not intended to be written to file, but instead
 * to be passed into the python environment.
 */

#ifndef __TRAJECTORYSINGLE_H_INCLUDED__
#define __TRAJECTORYSINGLE_H_INCLUDED__
#include "Observable.h"
#include "Particle.h"
#include "utils.h"
#include <vector>
#include <string>
#include <fstream>

class TrajectorySingle : public Observable
{
	public:
		/* shape: [time] [xyz]*/
		void record(
			std::vector<Particle> activeParticles,
			unsigned long int t);
		void writeBufferToFile();
};

#endif // __TRAJECTORYSINGLE_H_INCLUDED__
