/* ParticleTest.h
 * author: Christoph Froehner
 * 
 * Check if particle methods do what they should.
 */

#include <cxxtest/TestSuite.h>
#include "Particle.h"

class ParticleTest : public CxxTest::TestSuite
{
	public:
		void testMove(void)
		{
			Particle particle;
			std::array<double,3> x0 = {0.,0.,0.};
			std::array<double,3> dx = {1.,1.,1.};
			particle.position = x0;
			particle.move(dx);
			TS_ASSERT_EQUALS(particle.position[0],1.);
			TS_ASSERT_EQUALS(particle.position[1],1.);
			TS_ASSERT_EQUALS(particle.position[2],1.);
		}
};
