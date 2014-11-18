/* ActiveParticles.cpp
 * author: Christoph Froehner
 *
 */

#include "ActiveParticles.h"
#include <stdio.h>

// TODO this function inserts always at the beginning
void ActiveParticles::addParticle(Particle * particle)
{
	std::vector<Particle>::iterator it = this->container.begin();
	this->container.insert(it, *particle);
}
