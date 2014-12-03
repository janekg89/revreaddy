/* RandomTest.h
 * author: Christoph Froehner
 * 
 * Check if random methods do what they should.
 */

#include <cxxtest/TestSuite.h>
#include "Random.h"
#include <vector>
#include <string>

class RandomTest : public CxxTest::TestSuite
{
	public:
		/* construct 3 types of randomGenerators
		 * and check if they dont produce same
		 * numbers.
		 */
		void test_constructor_differentTypes(void)
		{
			Random randomRanlxs0("ranlxs0");
			Random randomTaus("taus");
			Random randomMt("mt19937");
			double x1 = randomRanlxs0.normal();
			double x2 = randomTaus.normal();
			double x3 = randomMt.normal();
			TS_ASSERT_DIFFERS(x1, x2);
			TS_ASSERT_DIFFERS(x2, x3); 
		}

		/* check if the three types of random generators
		 * produce normal distributed numbers that are
		 * smaller than 10 (abs. value). 
		 * Although there is a small chance to draw a large number
		 */
		void test_normal_inRange(void)
		{
			std::vector<std::string> types;
			types.push_back("ranlxs0");
			types.push_back("taus");
			types.push_back("mt19937");
			for (auto&& type : types)
			{
				Random * random = new Random(type);
				double x = random->normal();
				TS_ASSERT(x > -10.);
				TS_ASSERT(x <  10.);
				delete random;
			}
		}
};
