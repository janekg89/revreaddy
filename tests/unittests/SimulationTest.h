/* SimulationTest.h
 * 
 * A lot of unit tests of the Simulation class.
 */

#include <cxxtest/TestSuite.h>
// Simulation.h basically includes every other class used.
#include "Simulation.h"
#include "utils.h"

class SimulationTest : public CxxTest::TestSuite
{
	public:
		
		/* Add new particle types and check if this was applied
		 * correctly. */
		void test_newType(void)
		{
			Simulation * sim = new Simulation();
			// name, radius, diffConst, reactionRadius
			sim->new_Type("testtype", 2., 3., 4.);
			sim->new_Type("foo", 4., 5., 6.);
			// check the size of typeDict
			TS_ASSERT_EQUALS(sim->typeDict->names.size(), 2);
			TS_ASSERT_EQUALS(sim->typeDict->radii.size(), 2);
			TS_ASSERT_EQUALS(sim->typeDict->diffusionConstants.size(), 2);
			TS_ASSERT_EQUALS(sim->typeDict->reactionRadii.size(), 2);
			TS_ASSERT_EQUALS(sim->typeDict->getNumberOfTypes(), 2);

			// check for the values of typeDict
			TS_ASSERT_EQUALS(sim->typeDict->names[0], "testtype");
			TS_ASSERT_EQUALS(sim->typeDict->names[1], "foo");
			TS_ASSERT_DIFFERS(sim->typeDict->names[1], " foo");
			
			TS_ASSERT_EQUALS(sim->typeDict->radii[0], 2.);
			TS_ASSERT_EQUALS(sim->typeDict->radii[1], 4.);
			TS_ASSERT_DIFFERS(sim->typeDict->radii[1], 10.);

			TS_ASSERT_EQUALS(sim->typeDict->diffusionConstants[0], 3.);
			TS_ASSERT_EQUALS(sim->typeDict->diffusionConstants[1], 5.);
			TS_ASSERT_DIFFERS(sim->typeDict->diffusionConstants[1], 10.);

			TS_ASSERT_EQUALS(sim->typeDict->reactionRadii[0], 4.);
			TS_ASSERT_EQUALS(sim->typeDict->reactionRadii[1], 6.);
			TS_ASSERT_DIFFERS(sim->typeDict->reactionRadii[1], 10.);
			
			delete sim;
		}
		
		/* Add two particles and check whether their positions and
		 * types are applied correctly. */
		void test_addParticle(void)
		{
			Simulation * sim = new Simulation();
			sim->new_Type("soft1", 2., 3., 4.);
			sim->new_Type("soft2", 5., 6., 7.);
			std::vector<double> x1 = {1.,2.,3.};
			std::vector<double> x2 = {-3.,-2.,-1.};
			sim->addParticle(x1, 0);
			sim->addParticle(x2, 1);
			std::vector<double> x3 = {1.,2.,3.};
			std::vector<double> x4 = {-3.,-2.,-1.};
			TS_ASSERT_EQUALS(sim->activeParticles[0].position, x3);
			TS_ASSERT_EQUALS(sim->activeParticles[1].position, x4);
			TS_ASSERT_EQUALS(sim->activeParticles[0].typeId, 0);
			TS_ASSERT_EQUALS(sim->activeParticles[1].typeId, 1);
			delete sim;
		}
		
		/* Add a particle with diffusionConstant = 0. outside the box
		 * at {6.,0.,0.}. Boxsize is 10.
		 * After one propagate() it should be at {-4.,0.,0.}
		 */
		void test_propagate_periodicBoundaryCondition(void)
		{
			Simulation * sim       = new Simulation();
			sim->isPeriodic        = true;
			sim->boxsize           = 10.;
			sim->temperature       = 1.;
			sim->kBoltzmann        = 1.;
			std::vector<double> x0 = {6.,0.,0.};
			sim->new_Type("soft", 1., 0., 1.);
			sim->addParticle(x0, 0);
			sim->propagate();
			std::vector<double> x1 = {-4.,0.,0.};
			TS_ASSERT_EQUALS(sim->activeParticles[0].position, x1);
			delete sim;
		}

		/* Place two particles central in the box so that minimum
		 * distance is euclidean distance.
		 * Also place two particles on the edges of the box so that
		 * minimum distance is the distance "through the wall".
		 */
		void test_getMinDistance_minimumImageConvention(void)
		{
			Simulation * sim         = new Simulation();
			sim->isPeriodic          = true;
			sim->boxsize             = 10.;
			std::vector<double> x1  = {1.,1.,1.};
			std::vector<double> x2  = {-1.,-1.,-1.};
			std::vector<double> r12 = {-2.,-2.,-2.}; 
			std::vector<double> r   = {0.,0.,0.};
			getMinDistanceVector(
				r,
				x1,
				x2,
				sim->isPeriodic,
				sim->boxsize);
			TS_ASSERT_EQUALS(r, r12);
			std::vector<double> x3  = {4.,4.,4.};
			std::vector<double> x4  = {-4.,-4.,-4.};
			std::vector<double> r34 = {2.,2.,2.};
			getMinDistanceVector(
				r,
				x3,
				x4,
				sim->isPeriodic,
				sim->boxsize);
			TS_ASSERT_EQUALS(r, r34);
			delete sim;
		}
		
		/* Given 2 particles, check that the calculated force is correctly
		 * assigned to each particle. One particle is at {1,1,1}, the other
		 * at {2,3,4}. The force constant is k=1 and the particles' radii
		 * are 2. so that the cutoff is s=4. The force is
		 * F = 2* k * (1 - s/R) * r_12, where r_12 is the vector pointing from
		 * particle 1 to particle 2. The energy is E = k*(R - s)**2
		 */
		void test_calculateSingleForceEnergy_twoSoftParticles(void)
		{
			Simulation * sim = new Simulation();
			sim->isPeriodic = true;
			sim->boxsize = 10.;
			sim->new_Type("soft", 2., 0., 1.);
			sim->new_SoftRepulsion("softy", {0, 0}, 1.);
			std::vector<double> x0 = {1.,1.,1.};
			sim->addParticle(x0, 0);
			x0 = {2.,3.,4.};
			sim->addParticle(x0, 0);
			//printArray(sim->activeParticles[0].cumulativeForce);
			sim->calculateSingleForceEnergy(0, 1);
			//printArray(sim->activeParticles[0].cumulativeForce);
			std::vector<double> f = {1.,2.,3.};
			double preFactor = -0.13808993529939517;
			f[0] *= preFactor;
			f[1] *= preFactor;
			f[2] *= preFactor;
			TS_ASSERT_DELTA(sim->activeParticles[0].cumulativeForce, f, 0.00000001);
			f[0] *= -1.;
			f[1] *= -1.;
			f[2] *= -1.;
			TS_ASSERT_DELTA(sim->activeParticles[1].cumulativeForce, f, 0.00000001);
			TS_ASSERT_DELTA(sim->energy, 0.066740905808468948, 0.0000001);
			delete sim;
		}

		/* create a new interaction. check the cutoff distance. */
		void test_new_SoftRepulsionAndLennardJones(void)
		{
			Simulation * sim = new Simulation();
			sim->new_Type("A", 2., 1., 1.);
			sim->new_Type("B", 5., 1., 1.);
			sim->new_SoftRepulsion("rep", {0, 1}, 2.);
			TS_ASSERT_EQUALS(sim->possibleForces[0]->cutoff, 7.);
			sim->new_SoftRepulsion("rep2", {0, 0}, 3.);
			TS_ASSERT_EQUALS(sim->possibleForces[1]->cutoff, 4.);
			sim->new_LennardJones("lj", {1, 1}, 4.);
			TS_ASSERT_EQUALS(sim->possibleForces[2]->cutoff, 25.);
		}

		/* create two systems, initially the same. calculate forces
		 * in one system with neighborlattice and in the other
		 * without neighborlattice. Results (forces and energy)
		 * should be EXACTLY the same. */
		void test_neighborLattice_ForcesEnergiesSameAsWithout_soft(void)
		{
			Simulation * sim1 = new Simulation();
			Simulation * sim2 = new Simulation();
			sim1->boxsize = 10.;
			sim2->boxsize = 10.;
			sim1->new_Type("A", .5, 1., 1.);
			sim2->new_Type("A", .5, 1., 1.);
			sim1->new_SoftRepulsion("A<->A", {0,0}, 1.);
			sim2->new_SoftRepulsion("A<->A", {0,0}, 1.);
			sim1->addParticle({-4.,-3.,0.}, 0);
			sim2->addParticle({-4.,-3.,0.}, 0);
			sim1->addParticle({-3.3,-3.,0.}, 0);
			sim2->addParticle({-3.3,-3.,0.}, 0);
			sim1->addParticle({-1.9,-2.1,0.}, 0);
			sim2->addParticle({-1.9,-2.1,0.}, 0);
			sim1->addParticle({-2.2,-1.8,0.}, 0);
			sim2->addParticle({-2.2,-1.8,0.}, 0);
			
			sim1->calculateInteractionForcesEnergiesWithLattice(10);
			sim2->calculateInteractionForcesEnergiesNaive();

			TS_ASSERT_EQUALS(sim1->energy, sim2->energy);
			for (unsigned int i=0; i<4; i++) {
				TS_ASSERT_EQUALS(
					sim1->activeParticles[i].cumulativeForce,
					sim2->activeParticles[i].cumulativeForce);
			}
		}

		/* same as before but with LennardJones particles */
		void test_neighborLattice_ForcesEnergiesSameAsWithout_lj(void)
		{
			Simulation * sim1 = new Simulation();
			Simulation * sim2 = new Simulation();
			sim1->boxsize = 10.;
			sim2->boxsize = 10.;
			sim1->new_Type("A", .5, 1., 1.);
			sim2->new_Type("A", .5, 1., 1.);
			sim1->new_LennardJones("A<->A", {0,0}, 1.);
			sim2->new_LennardJones("A<->A", {0,0}, 1.);
			sim1->addParticle({-4.2,-3.,0.}, 0);
			sim2->addParticle({-4.2,-3.,0.}, 0);
			sim1->addParticle({-3.3,-3.,0.}, 0);
			sim2->addParticle({-3.3,-3.,0.}, 0);
			sim1->addParticle({-1.8,-2.1,0.}, 0);
			sim2->addParticle({-1.8,-2.1,0.}, 0);
			sim1->addParticle({-2.5,-1.6,0.}, 0);
			sim2->addParticle({-2.5,-1.6,0.}, 0);
			
			sim1->calculateInteractionForcesEnergiesWithLattice(4);
			sim2->calculateInteractionForcesEnergiesNaive();

			TS_ASSERT_EQUALS(sim1->energy, sim2->energy);
			for (unsigned int i=0; i<4; i++) {
				TS_ASSERT_EQUALS(
					sim1->activeParticles[i].cumulativeForce,
					sim2->activeParticles[i].cumulativeForce);
			}
		}

};
