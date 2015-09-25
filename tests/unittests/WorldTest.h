/* WorldTest.h */

#include <cxxtest/TestSuite.h>
#include "World.h"

class WorldTest : public CxxTest::TestSuite
{
	public:

	/* Add two particles and check whether their positions and
	 * types are applied correctly. */
	void test_addParticle(void)
	{
		World * world = new World();
		std::vector<double> x1 = {1.,2.,3.};
		std::vector<double> x2 = {-3.,-2.,-1.};
		world->addParticle(x1, 0);
		world->addParticle(x2, 1);
		std::vector<double> x3 = {1.,2.,3.};
		std::vector<double> x4 = {-3.,-2.,-1.};
		TS_ASSERT_EQUALS(world->activeParticles[0].position, x3);
		TS_ASSERT_EQUALS(world->activeParticles[1].position, x4);
		TS_ASSERT_EQUALS(world->activeParticles[0].typeId, 0);
		TS_ASSERT_EQUALS(world->activeParticles[1].typeId, 1);
		delete world;
	}	
};