/* ActiveParticles.h
 * author: Christoph Froehner
 *
 * This is the container that will hold each particle.
 * It will derive also 'consideredParticles', which
 * will hold all particles that are actually propagated.
 *
 */

#ifndef __ACTIVEPARTICLES_H_INCLUDED__
#define __ACTIVEPARTICLES_H_INCLUDED__
#include <vector>
#include "Particle.h"

class ActiveParticles
{
	//private:

	public:
		std::vector<Particle> container;
		void calculateUnaryForces();
		// idea: loop over particlepair (i,j): 0<i<n and i<j<n
		// lookup particletype: determines interaction
		// also lookup group of particle
		// question: should grouped particles interact with exclusive
		// forces and standard potentials with all others?
		// maybe extra function: calculateGroupForces();
		void calculateBinaryForces();

		// will be needed for adaptive timestepping
		void constructConsideredParticles();

		void addParticle(Particle * particle);

};

#endif
