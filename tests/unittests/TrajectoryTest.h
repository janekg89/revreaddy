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
			Trajectory traj("test.txt");
			std::vector<Particle> activeParticles;
			std::vector<double> x0 = {1.,2.,3.};
			Particle p1;
			Particle p2;
			p1.position = x0;
			p2.position = x0;
			p1.typeId = 0;
			p2.typeId = 1;
			activeParticles.push_back(p1);
			activeParticles.push_back(p2);
			activeParticles.push_back(p1);
			double t = 42.;
			traj.record(activeParticles, t);
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
