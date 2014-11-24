/* SingleParticleDiffusion.cpp
 * author: Christoph Froehner
 */

#include "SingleParticleDiffusion.h"
#include "utils.h"

void SingleParticleDiffusion::recordObservables(unsigned long int t)
{
	squaredDistances[t] = squaredDistance(initialPosition, activeParticles[0].position);
	trajectory[t] = activeParticles[0].position;
}
