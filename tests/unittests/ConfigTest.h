/* ConfigTest.h */

#include <cxxtest/TestSuite.h>
#include "Config.h"

class ConfigTest : public CxxTest::TestSuite
{
	public:

	/* Add new particle types and check if this was applied
	 * correctly. */
	void test_newType(void)
	{
		World * world = new World();
		Config * con = new Config(world);
		// name, radius, diffConst, reactionRadius
		con->new_Type("testtype", 2., 3., 4.);
		con->new_Type("foo", 4., 5., 6.);
		// check the size of typeDict
		TS_ASSERT_EQUALS(con->typeDict.size(), 2);

		// check for the values of typeDict
		TS_ASSERT_EQUALS(con->typeDict[0].name, "testtype");
		TS_ASSERT_EQUALS(con->typeDict[1].name, "foo");
		TS_ASSERT_DIFFERS(con->typeDict[1].name, " foo");
		
		TS_ASSERT_EQUALS(con->typeDict[0].radius, 2.);
		TS_ASSERT_EQUALS(con->typeDict[1].radius, 4.);
		TS_ASSERT_DIFFERS(con->typeDict[1].radius, 10.);

		TS_ASSERT_EQUALS(con->typeDict[0].diffusionConstant, 3.);
		TS_ASSERT_EQUALS(con->typeDict[1].diffusionConstant, 5.);
		TS_ASSERT_DIFFERS(con->typeDict[1].diffusionConstant, 10.);

		TS_ASSERT_EQUALS(con->typeDict[0].reactionRadius, 4.);
		TS_ASSERT_EQUALS(con->typeDict[1].reactionRadius, 6.);
		TS_ASSERT_DIFFERS(con->typeDict[1].reactionRadius, 10.);
		
		delete con;
		delete world;
	}

	/* create a new interaction. check the cutoff distance. */
	void test_new_SoftRepulsionAndLennardJones(void)
	{
		World * world = new World();
		Config * config = new Config(world);
		config->new_Type("A", 2., 1., 1.);
		config->new_Type("B", 5., 1., 1.);
		config->new_SoftRepulsion("rep", {0, 1}, 2.);
		TS_ASSERT_EQUALS(config->possibleInteractions[0]->cutoff, 7.);
		config->new_SoftRepulsion("rep2", {0, 0}, 3.);
		TS_ASSERT_EQUALS(config->possibleInteractions[1]->cutoff, 4.);
		config->new_LennardJones("lj", {1, 1}, 4.);
		TS_ASSERT_EQUALS(config->possibleInteractions[2]->cutoff, 25.);
		delete config;
		delete world;
	}
};