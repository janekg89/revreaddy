/* ActiveParticles.cpp
 * author: Christoph Froehner
 *
 */

#include "ActiveParticles.h"
#include <stdio.h>

// TODO this function inserts always at the beginning
void ActiveParticles::addParticle(Particle * particle)
{
	this->container.push_back(*particle);
}
