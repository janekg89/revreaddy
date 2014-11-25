/* GenericSimulation.h
 * author: Christoph Froehner
 * 
 * Provide a generic class for performing all
 * sorts of particle simulation.
 */

#ifndef __GENERICSIMULATION_H_INCLUDED__
#define __GENERICSIMULATION_H_INCLUDED__
#include "Simulation.h"

// TODO think of a way to record whatever observable the user
// wants to record, from reading the 'observables' file.
class GenericSimulation: public Simulation
{
	public:
		void recordObservables(unsigned long int t);
};
#endif
