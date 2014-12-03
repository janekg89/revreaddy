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
		void test_move_usualOperation(void)
		{
			Particle particle;
			std::array<double,3> x0 = {0.,0.,0.};
			std::array<double,3> dx = {1.,1.,1.};
			std::array<double,3> x1 = {1.,1.,1.};
			particle.position = x0;

			particle.move(dx);

			TS_ASSERT_EQUALS(particle.position, x1);
		}

		void test_addForce_usualOperation(void)
		{
			Particle particle;
			std::array<double,3> force = {1.,-2.,3.};
			std::array<double,3> res = {1.,-2.,3.};

			particle.addForce(force);

			TS_ASSERT_EQUALS(particle.cumulativeForce, res);
		}

		void test_resetForce_usualOperation(void)
		{
			Particle particle;
			std::array<double,3> force = {2.,1.,-3.};
			std::array<double,3> zero = {0.,0.,0.};
			particle.cumulativeForce = force;

			particle.resetForce();

			TS_ASSERT_EQUALS(particle.cumulativeForce, zero);
		}
};
