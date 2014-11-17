/* SingleParticleDiffusion.cpp
 * author: Christoph Froehner
 */

#include "SingleParticleDiffusion.h"
#include "utils.h"

void SingleParticleDiffusion::run(ActiveParticles * activeParticles, Random * random)
{
	double initialPosition[3] = {0., 0., 0.};
	activeParticles->container[0].position[0] = initialPosition[0];
	activeParticles->container[0].position[1] = initialPosition[1];
	activeParticles->container[0].position[2] = initialPosition[2];
	this->meanSquaredDistances[0] = 0.;
	for (int t = 1; t < this->maxTime; t++)
	{
		this->propagate(activeParticles, random);
		this->meanSquaredDistances[t] = squaredDistance(initialPosition, activeParticles->container[0].position);
	}
}
