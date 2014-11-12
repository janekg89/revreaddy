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

class ActiveParticles
{
	private:
		std::vector<Particle> container;
};

#endif
