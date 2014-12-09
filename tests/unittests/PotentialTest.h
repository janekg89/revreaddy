/* PotentialTest.cpp */

#include <cxxtest/TestSuite.h>
#include "Potential.h"

class PotentialTest : public CxxTest::TestSuite
{
	public:
		/***
		 * One particle at {1,1,1} and another at {2,3,4}, so that
		 * r_ij = {1,2,3}. The force is F = k * (1 - s/R) * r_ij,
		 * where k is 2, R is sqrt(14) and s is the cutoff (sum of
		 * the particles' radii).
		 * The force should only be repulsive, so that 
		 * F = {0,0,0} if R > s.
		 ***/
		void test_softcoreForce_usualOperation(void)
		{
			Potential potential;
			std::array<double,3> r_ij = {1.,2.,3.};
			std::array<double,3> force;
			force = potential.softcoreForce(r_ij, 14., 16., 2.);
			std::array<double,3> expForce;
			double expPrefactor = -0.13808993529939517;
			expForce[0] = expPrefactor * r_ij[0];
			expForce[1] = expPrefactor * r_ij[1];
			expForce[2] = expPrefactor * r_ij[2];
			TS_ASSERT_EQUALS(force, expForce);
			r_ij = {2.,4.,5.};
			force = potential.softcoreForce(r_ij, 45., 16., 2.);
			expForce = {0.,0.,0.};
			TS_ASSERT_EQUALS(force, expForce);
		}
};
