/* SingleParticleDiffusion.cpp
 * author: Christoph Froehner
 */

#include "SingleParticleDiffusion.h"
#include "utils.h"

void SingleParticleDiffusion::run(ActiveParticles * activeParticles, Random * random)
{
	std::array<double, 3> initialPosition = {0., 0., 0.};
	activeParticles->container[0].position = initialPosition;
	this->squaredDistances[0] = 0.;
	for (int t = 1; t < this->maxTime; t++)
	{
		this->propagate(activeParticles, random);
		this->squaredDistances[t] = squaredDistance(initialPosition, activeParticles->container[0].position);
	}
}
