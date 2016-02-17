/* TrajectoryTest.h
 * author: Christoph Froehner
 * 
 * Check if Observable methods do what they should.
 */

#include <cxxtest/TestSuite.h>
// Trajectory.h is a child of Observable and includes Particle.h
#include "Trajectory.h"

class TrajectoryTest : public CxxTest::TestSuite
{
	public:
		void test_record_usualOperation(void)
		{
			World * world = new World();
			Trajectory traj(1,1,"test.txt");
			std::vector<double> x0 = {1.,2.,3.};
			world->addParticle(x0, 0);
			world->addParticle(x0, 1);
			world->addParticle(x0, 0);
			double t = 42.;
			traj.record(world, t);
			// check if only one timestep is in the trajectory
			TS_ASSERT_EQUALS(traj.trajectory.size(), 1);
			// check if there are 3 particles in the first trajectory entry
			TS_ASSERT_EQUALS(traj.trajectory[0].size(), 3);
			// check if the particles' time is indeed 42
			TS_ASSERT_EQUALS(traj.trajectory[0][0].particleTime, 42.);
			TS_ASSERT_EQUALS(traj.trajectory[0][1].particleTime, 42.);
			TS_ASSERT_EQUALS(traj.trajectory[0][2].particleTime, 42.);
			// check if the particles' position is {1.,2.,3.}
			std::vector<double> x1 = {1.,2.,3.};
			TS_ASSERT_EQUALS(traj.trajectory[0][0].particleCoordinates, x1);
			TS_ASSERT_EQUALS(traj.trajectory[0][1].particleCoordinates, x1);
			TS_ASSERT_EQUALS(traj.trajectory[0][2].particleCoordinates, x1);
			// check if the particles' type is 0 or 1
			TS_ASSERT_EQUALS(traj.trajectory[0][0].particleTypeId, 0);
			TS_ASSERT_EQUALS(traj.trajectory[0][1].particleTypeId, 1);
			TS_ASSERT_EQUALS(traj.trajectory[0][2].particleTypeId, 0);
		}
};
