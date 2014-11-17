/* ActiveParticles.cpp
 * author: Christoph Froehner
 *
 */

#include "ActiveParticles.h"
#include <stdio.h>

void ActiveParticles::addParticle(Particle * particle)
{
	std::vector<Particle>::iterator it = this->container.begin();
	this->container.insert(it, *particle);
}
