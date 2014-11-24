/* SingleParticleDiffusion.h
 * author: Christoph Froehner
 *
 * The aim of this simulation is to prove obediance
 * to the diffusion equation.
 */

#ifndef __SINGLEPARTICLEDIFFUSION_H_INCLUDED__
#define __SINGLEPARTICLEDIFFUSION_H_INCLUDED__
#include "Simulation.h"
#include <array>
#include <vector>

class SingleParticleDiffusion: public Simulation
{
	public:
		std::vector<double> squaredDistances;
		std::array<double, 3> initialPosition;
		void recordObservables(unsigned long int t);
		std::vector< std::array<double, 3> > trajectory;
};

#endif
