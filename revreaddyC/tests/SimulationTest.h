/* TODO googletest */
/* SimulationTest.h */
#include "Simulation.h"
#include "utils.h"

class SimulationTest : public CxxTest::TestSuite
{
	public:
	
	/* Add a particle with diffusionConstant = 0. outside the box
	 * at {6.,0.,0.}. Boxsize is 10.
	 * After one propagate() it should be at {-4.,0.,0.}
	 */
	void test_propagate_periodicBoundaryCondition(void)
	{
		Simulation * sim = new Simulation();
		World * world = sim->world;
		Config * config = sim->config;
		config->isPeriodic = true;
		config->boxsize           = 10.;
		config->temperature       = 1.;
		config->kBoltzmann        = 1.;
		std::vector<double> x0 = {6.,0.,0.};
		config->new_Type("soft", 1., 0., 1.);
		world->addParticle(x0, 0);
		sim->propagateDynamics();
		std::vector<double> x1 = {-4.,0.,0.};
		TS_ASSERT_EQUALS(world->activeParticles[0].position, x1);
		delete sim;
	}

	/* Place two particles central in the box so that minimum
	 * distance is euclidean distance.
	 * Also place two particles on the edges of the box so that
	 * minimum distance is the distance "through the wall".
	 */
	void test_getMinDistance_minimumImageConvention(void)
	{
		Simulation * sim        = new Simulation();
		sim->config->isPeriodic = true;
		sim->config->boxsize    = 10.;
		std::vector<double> x1  = {1.,1.,1.};
		std::vector<double> x2  = {-1.,-1.,-1.};
		std::vector<double> r12 = {-2.,-2.,-2.}; 
		std::vector<double> r   = {0.,0.,0.};
		getMinDistanceVector(
			r,
			x1,
			x2,
			sim->config->isPeriodic,
			sim->config->boxsize);
		TS_ASSERT_EQUALS(r, r12);
		std::vector<double> x3  = {4.,4.,4.};
		std::vector<double> x4  = {-4.,-4.,-4.};
		std::vector<double> r34 = {2.,2.,2.};
		getMinDistanceVector(
			r,
			x3,
			x4,
			sim->config->isPeriodic,
			sim->config->boxsize);
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
		sim->config->isPeriodic = true;
		sim->config->boxsize = 10.;
		sim->config->new_Type("soft", 2., 0., 1.);
		sim->config->new_SoftRepulsion("softy", {0, 0}, 1.);
		std::vector<double> x0 = {1.,1.,1.};
		sim->world->addParticle(x0, 0);
		x0 = {2.,3.,4.};
		sim->world->addParticle(x0, 0);
		//printArray(sim->activeParticles[0].cumulativeForce);
		sim->calculateSingleForceEnergy(0, 1);
		//printArray(sim->activeParticles[0].cumulativeForce);
		std::vector<double> f = {1.,2.,3.};
		double preFactor = -0.13808993529939517;
		f[0] *= preFactor;
		f[1] *= preFactor;
		f[2] *= preFactor;
		TS_ASSERT_DELTA(sim->world->activeParticles[0].cumulativeForce, f, 1e-9);
		f[0] *= -1.;
		f[1] *= -1.;
		f[2] *= -1.;
		TS_ASSERT_DELTA(sim->world->activeParticles[1].cumulativeForce, f, 1e-9);
		TS_ASSERT_DELTA(sim->world->energy, 0.066740905808468948, 1e-9);
		delete sim;
	}

	/* Two types of particles A and B with reactionRadii 2 and 1
	 * respectively. After initial setup the activePairs should
	 * contain the correct pairs of indices */
	void test_activePairs_correctEntries(void)
	{
		Simulation sim;
		sim.config->new_Type("A", 1., 1., 2.);
		sim.config->new_Type("B", 1., 1., 1.);
		sim.config->boxsize = 10.;
		sim.world->addParticle({0., -2., 0.}, 0); // Type A uniqueId 0
		sim.world->addParticle({2., -2., 0.}, 1); // B 1
		sim.world->addParticle({0., 3., 0.}, 1);  // B 2
		sim.world->addParticle({2., 3., 0.}, 0);  // A 3
		sim.world->addParticle({3., 3., 0.}, 1);  // B 4
		sim.calculateInteractionForcesEnergiesNaive();
		TS_ASSERT_EQUALS(sim.world->activePairs.size(), 3);
		std::vector<unsigned long int> pair1 = {0, 1};
		std::vector<unsigned long int> pair2 = {2, 3};
		std::vector<unsigned long int> pair3 = {3, 4};
		TS_ASSERT_EQUALS(sim.world->activePairs[0], pair1);
		TS_ASSERT_EQUALS(sim.world->activePairs[1], pair2);
		TS_ASSERT_EQUALS(sim.world->activePairs[2], pair3);
	}

	/* Create a current state (energy, activeParticles, activePairs)
	 * and check if it is saved correctly by comparing individual 
	 * elements of current and old variables*/
	void test_saveOldState(void)
	{
		Simulation sim;
		sim.config->new_Type("A", 1., 1., 2.);
		sim.config->new_Type("B", 1., 1., 1.);
		sim.config->boxsize = 10.;
		/* Set current state */
		sim.world->addParticle({0., -2., 0.}, 0); // Type A uniqueId 0
		sim.world->addParticle({2., -2., 0.}, 1); // B 1
		sim.world->addParticle({0., 3., 0.}, 1);  // B 2
		sim.world->addParticle({2., 3., 0.}, 0);  // A 3
		sim.world->addParticle({3., 3., 0.}, 1);  // B 4
		std::vector<unsigned long int> pair1 = {0, 1};
		std::vector<unsigned long int> pair2 = {2, 3};
		std::vector<unsigned long int> pair3 = {3, 4};
		sim.world->activePairs.push_back(pair1);
		sim.world->activePairs.push_back(pair2);
		sim.world->activePairs.push_back(pair3);
		sim.world->energy = 42.;
		/* Perform save */
		sim.saveOldState();
		/* Check for equality of old and current variables.
		 * First check activeParticles. */
		TS_ASSERT_EQUALS(
			sim.world->oldActiveParticles.size(), sim.world->activeParticles.size());
		for (unsigned int i=0; i<sim.world->activeParticles.size(); i++) {
			TS_ASSERT_EQUALS(
				sim.world->oldActiveParticles[i], sim.world->activeParticles[i]);
		}
		
		TS_ASSERT_EQUALS(
			sim.world->oldActivePairs.size(), sim.world->activePairs.size());
		TS_ASSERT_EQUALS(sim.world->oldActivePairs[0], pair1);
		TS_ASSERT_EQUALS(sim.world->oldActivePairs[1], pair2);
		TS_ASSERT_EQUALS(sim.world->oldActivePairs[2], pair3);
		TS_ASSERT_EQUALS(sim.world->oldEnergy, 42.);
	}

	/* Same as test_saveOldState but testing restoreOldState */
	void test_restoreOldState(void)
	{
		Simulation sim;
		sim.config->new_Type("A", 1., 1., 2.);
		sim.config->new_Type("B", 1., 1., 1.);
		sim.config->boxsize = 10.;
		/* Set current state */
		sim.world->addParticle({0., -2., 0.}, 0); // Type A uniqueId 0
		sim.world->addParticle({2., -2., 0.}, 1); // B 1
		sim.world->addParticle({0., 3., 0.}, 1);  // B 2
		sim.world->addParticle({2., 3., 0.}, 0);  // A 3
		sim.world->addParticle({3., 3., 0.}, 1);  // B 4
		sim.world->oldActiveParticles = sim.world->activeParticles;
		sim.world->activeParticles.clear();
		std::vector<unsigned long int> pair1 = {0, 1};
		std::vector<unsigned long int> pair2 = {2, 3};
		std::vector<unsigned long int> pair3 = {3, 4};
		sim.world->oldActivePairs.push_back(pair1);
		sim.world->oldActivePairs.push_back(pair2);
		sim.world->oldActivePairs.push_back(pair3);
		sim.world->oldEnergy = 42.;
		/* Perform restore */
		sim.restoreOldState();
		/* Check for equality of old and current variables.
		 * First check activeParticles. */
		TS_ASSERT_EQUALS(
			sim.world->oldActiveParticles.size(), sim.world->activeParticles.size());
		for (unsigned int i=0; i<sim.world->activeParticles.size(); i++) {
			TS_ASSERT_EQUALS(
				sim.world->oldActiveParticles[i], sim.world->activeParticles[i]);
		}
		
		TS_ASSERT_EQUALS(
			sim.world->oldActivePairs.size(), sim.world->activePairs.size());
		TS_ASSERT_EQUALS(sim.world->activePairs[0], pair1);
		TS_ASSERT_EQUALS(sim.world->activePairs[1], pair2);
		TS_ASSERT_EQUALS(sim.world->activePairs[2], pair3);
		TS_ASSERT_EQUALS(sim.world->energy, 42.);		
	}

	/* create two systems, initially the same. calculate forces
	 * in one system with neighborlattice and in the other
	 * without neighborlattice. Results (forces and energy)
	 * should be the same. */
	void test_neighborLattice_ForcesEnergiesSameAsWithout_soft(void)
	{
		Simulation * sim1 = new Simulation();
		Simulation * sim2 = new Simulation();
		sim1->config->boxsize = 10.;
		sim2->config->boxsize = 10.;
		sim1->config->new_Type("A", 1.1, 1., 1.);
		sim2->config->new_Type("A", 1.1, 1., 1.);
		sim1->config->new_SoftRepulsion("A<->A", {0,0}, 10.);
		sim2->config->new_SoftRepulsion("A<->A", {0,0}, 10.);
		
		sim1->world->addParticle({ 0.558377928961 , 0.397901875093 , -2.33439151343 },0);
		sim2->world->addParticle({ 0.558377928961 , 0.397901875093 , -2.33439151343 },0);
		sim1->world->addParticle({ 2.03345233286 , -4.84908651104 , 3.19673596552 },0);
		sim2->world->addParticle({ 2.03345233286 , -4.84908651104 , 3.19673596552 },0);
		sim1->world->addParticle({ 4.27250843873 , -4.13179374924 , 0.113051639401 },0);
		sim2->world->addParticle({ 4.27250843873 , -4.13179374924 , 0.113051639401 },0);
		sim1->world->addParticle({ 2.78244058766 , -2.21138694369 , -2.13996550751 },0);
		sim2->world->addParticle({ 2.78244058766 , -2.21138694369 , -2.13996550751 },0);
		sim1->world->addParticle({ 3.34142436673 , 2.71053025753 , -1.83795532374 },0);
		sim2->world->addParticle({ 3.34142436673 , 2.71053025753 , -1.83795532374 },0);
		sim1->world->addParticle({ 3.17912954464 , -2.39969579008 , 1.90934936849 },0);
		sim2->world->addParticle({ 3.17912954464 , -2.39969579008 , 1.90934936849 },0);
		sim1->world->addParticle({ -2.0858953566 , -3.24876196998 , 4.40617622431 },0);
		sim2->world->addParticle({ -2.0858953566 , -3.24876196998 , 4.40617622431 },0);
		sim1->world->addParticle({ -2.7347226244 , -3.59318561718 , 1.89550416966 },0);
		sim2->world->addParticle({ -2.7347226244 , -3.59318561718 , 1.89550416966 },0);
		sim1->world->addParticle({ 3.40016731125 , 3.19962010176 , 3.73811860165 },0);
		sim2->world->addParticle({ 3.40016731125 , 3.19962010176 , 3.73811860165 },0);
		sim1->world->addParticle({ -3.66264291085 , -1.99237375949 , -4.25245575687 },0);
		sim2->world->addParticle({ -3.66264291085 , -1.99237375949 , -4.25245575687 },0);
		sim1->world->addParticle({ 0.281748081768 , 2.56389859258 , -2.87479159471 },0);
		sim2->world->addParticle({ 0.281748081768 , 2.56389859258 , -2.87479159471 },0);
		sim1->world->addParticle({ 0.318233298416 , -1.56094046143 , 2.28550044143 },0);
		sim2->world->addParticle({ 0.318233298416 , -1.56094046143 , 2.28550044143 },0);
		sim1->world->addParticle({ 1.29205921481 , 1.77426373706 , 3.55881762954 },0);
		sim2->world->addParticle({ 1.29205921481 , 1.77426373706 , 3.55881762954 },0);
		sim1->world->addParticle({ -3.67842677619 , 2.12158399814 , 4.81634458246 },0);
		sim2->world->addParticle({ -3.67842677619 , 2.12158399814 , 4.81634458246 },0);
		sim1->world->addParticle({ -4.23748061645 , -1.76192862683 , 2.87383521151 },0);
		sim2->world->addParticle({ -4.23748061645 , -1.76192862683 , 2.87383521151 },0);
		sim1->world->addParticle({ 4.66654217592 , 4.27524213312 , -4.43022193803 },0);
		sim2->world->addParticle({ 4.66654217592 , 4.27524213312 , -4.43022193803 },0);
		sim1->world->addParticle({ -3.83615941064 , -0.186641805865 , -2.10930566258 },0);
		sim2->world->addParticle({ -3.83615941064 , -0.186641805865 , -2.10930566258 },0);
		sim1->world->addParticle({ -4.29910730187 , 0.561186478851 , 3.19724755656 },0);
		sim2->world->addParticle({ -4.29910730187 , 0.561186478851 , 3.19724755656 },0);
		sim1->world->addParticle({ -1.04258264872 , 4.40034516842 , -4.20633487733 },0);
		sim2->world->addParticle({ -1.04258264872 , 4.40034516842 , -4.20633487733 },0);
		sim1->world->addParticle({ -0.483248084695 , 3.82190881754 , -1.16057823317 },0);
		sim2->world->addParticle({ -0.483248084695 , 3.82190881754 , -1.16057823317 },0);
		sim1->world->addParticle({ 2.77718015879 , -1.43586958375 , -0.654647659458 },0);
		sim2->world->addParticle({ 2.77718015879 , -1.43586958375 , -0.654647659458 },0);
		sim1->world->addParticle({ 3.87564798234 , -4.84347088855 , -2.1800763112 },0);
		sim2->world->addParticle({ 3.87564798234 , -4.84347088855 , -2.1800763112 },0);
		sim1->world->addParticle({ -0.76000036765 , 4.6861535406 , 2.93836596812 },0);
		sim2->world->addParticle({ -0.76000036765 , 4.6861535406 , 2.93836596812 },0);
		sim1->world->addParticle({ 0.926086632477 , -3.17256305255 , -0.551751586194 },0);
		sim2->world->addParticle({ 0.926086632477 , -3.17256305255 , -0.551751586194 },0);
		sim1->world->addParticle({ -0.466196801032 , 0.318874196587 , 1.25325756971 },0);
		sim2->world->addParticle({ -0.466196801032 , 0.318874196587 , 1.25325756971 },0);
		sim1->world->addParticle({ -1.9887082264 , -1.4143061935 , 2.42690296955 },0);
		sim2->world->addParticle({ -1.9887082264 , -1.4143061935 , 2.42690296955 },0);
		sim1->world->addParticle({ -2.66956475055 , 4.22390674538 , 0.494962078803 },0);
		sim2->world->addParticle({ -2.66956475055 , 4.22390674538 , 0.494962078803 },0);
		sim1->world->addParticle({ -1.37101029195 , -0.765030066283 , -4.56968767568 },0);
		sim2->world->addParticle({ -1.37101029195 , -0.765030066283 , -4.56968767568 },0);
		sim1->world->addParticle({ -2.05803582894 , 4.68885951471 , -2.06372233177 },0);
		sim2->world->addParticle({ -2.05803582894 , 4.68885951471 , -2.06372233177 },0);
		sim1->world->addParticle({ 1.64103181155 , -0.477802664316 , 3.43915170379 },0);
		sim2->world->addParticle({ 1.64103181155 , -0.477802664316 , 3.43915170379 },0);
		sim1->world->addParticle({ 1.40081742972 , 3.6420165099 , 0.363373031527 },0);
		sim2->world->addParticle({ 1.40081742972 , 3.6420165099 , 0.363373031527 },0);
		sim1->world->addParticle({ 1.89750867387 , 4.14033795384 , -3.70099875548 },0);
		sim2->world->addParticle({ 1.89750867387 , 4.14033795384 , -3.70099875548 },0);
		sim1->world->addParticle({ -3.96557621349 , -4.76987512754 , 3.39836665283 },0);
		sim2->world->addParticle({ -3.96557621349 , -4.76987512754 , 3.39836665283 },0);
		sim1->world->addParticle({ 0.45365214807 , 0.0785622878332 , -4.9810122132 },0);
		sim2->world->addParticle({ 0.45365214807 , 0.0785622878332 , -4.9810122132 },0);
		sim1->world->addParticle({ 2.90447763147 , 2.26559830081 , -3.86763647465 },0);
		sim2->world->addParticle({ 2.90447763147 , 2.26559830081 , -3.86763647465 },0);
		sim1->world->addParticle({ 0.553193655495 , -3.73911272994 , 2.0657931251 },0);
		sim2->world->addParticle({ 0.553193655495 , -3.73911272994 , 2.0657931251 },0);
		sim1->world->addParticle({ 1.25244674888 , -3.63964335237 , 4.94191835974 },0);
		sim2->world->addParticle({ 1.25244674888 , -3.63964335237 , 4.94191835974 },0);
		sim1->world->addParticle({ 0.358106318553 , -2.08280143372 , -3.91412153166 },0);
		sim2->world->addParticle({ 0.358106318553 , -2.08280143372 , -3.91412153166 },0);
		sim1->world->addParticle({ 4.79763033132 , 3.4219921353 , 2.09992572592 },0);
		sim2->world->addParticle({ 4.79763033132 , 3.4219921353 , 2.09992572592 },0);
		sim1->world->addParticle({ -4.57940646632 , 2.2485133573 , -2.99819767087 },0);
		sim2->world->addParticle({ -4.57940646632 , 2.2485133573 , -2.99819767087 },0);
		sim1->world->addParticle({ -2.23376267987 , 1.33378352796 , -4.43405751144 },0);
		sim2->world->addParticle({ -2.23376267987 , 1.33378352796 , -4.43405751144 },0);
		sim1->world->addParticle({ 2.04031156085 , -4.03804929517 , -2.9845386846 },0);
		sim2->world->addParticle({ 2.04031156085 , -4.03804929517 , -2.9845386846 },0);
		sim1->world->addParticle({ -2.83072035674 , 1.55165731998 , -0.732559276429 },0);
		sim2->world->addParticle({ -2.83072035674 , 1.55165731998 , -0.732559276429 },0);
		sim1->world->addParticle({ 1.93782561627 , -0.152819404501 , -3.82454080237 },0);
		sim2->world->addParticle({ 1.93782561627 , -0.152819404501 , -3.82454080237 },0);
		sim1->world->addParticle({ -1.89023186179 , 0.402099856928 , 3.50906551752 },0);
		sim2->world->addParticle({ -1.89023186179 , 0.402099856928 , 3.50906551752 },0);
		sim1->world->addParticle({ 3.81749035314 , 1.24581161939 , 4.60512596267 },0);
		sim2->world->addParticle({ 3.81749035314 , 1.24581161939 , 4.60512596267 },0);
		sim1->world->addParticle({ 0.490208816086 , 3.24820445642 , 4.619249877 },0);
		sim2->world->addParticle({ 0.490208816086 , 3.24820445642 , 4.619249877 },0);
		sim1->world->addParticle({ -2.1101617628 , 3.2918297016 , 4.31726696137 },0);
		sim2->world->addParticle({ -2.1101617628 , 3.2918297016 , 4.31726696137 },0);
		sim1->world->addParticle({ 4.50623223251 , -1.5994350353 , -2.53300859903 },0);
		sim2->world->addParticle({ 4.50623223251 , -1.5994350353 , -2.53300859903 },0);
		sim1->world->addParticle({ 1.43384096489 , 3.94040824581 , -1.82555768869 },0);
		sim2->world->addParticle({ 1.43384096489 , 3.94040824581 , -1.82555768869 },0);
		sim1->world->addParticle({ -3.45791511285 , 3.77490822303 , -1.36450109769 },0);
		sim2->world->addParticle({ -3.45791511285 , 3.77490822303 , -1.36450109769 },0);
		sim1->world->addParticle({ -3.31048693832 , -2.90788776886 , -0.485027921749 },0);
		sim2->world->addParticle({ -3.31048693832 , -2.90788776886 , -0.485027921749 },0);
		sim1->world->addParticle({ -2.26367401796 , -2.64904665948 , -3.02102158557 },0);
		sim2->world->addParticle({ -2.26367401796 , -2.64904665948 , -3.02102158557 },0);
		sim1->world->addParticle({ -4.06041900596 , -4.29136758251 , 1.32304654313 },0);
		sim2->world->addParticle({ -4.06041900596 , -4.29136758251 , 1.32304654313 },0);
		sim1->world->addParticle({ 3.63532594127 , -3.56456842938 , -4.83983885732 },0);
		sim2->world->addParticle({ 3.63532594127 , -3.56456842938 , -4.83983885732 },0);
		sim1->world->addParticle({ 4.26095358097 , -3.91328469304 , 3.15542985197 },0);
		sim2->world->addParticle({ 4.26095358097 , -3.91328469304 , 3.15542985197 },0);
		sim1->world->addParticle({ 0.641328733557 , -1.53677506212 , -2.34503393881 },0);
		sim2->world->addParticle({ 0.641328733557 , -1.53677506212 , -2.34503393881 },0);
		sim1->world->addParticle({ 2.84727254866 , 0.735355753206 , 0.696425378132 },0);
		sim2->world->addParticle({ 2.84727254866 , 0.735355753206 , 0.696425378132 },0);
		sim1->world->addParticle({ -3.88709759642 , 0.194476864104 , 0.782265067034 },0);
		sim2->world->addParticle({ -3.88709759642 , 0.194476864104 , 0.782265067034 },0);
		sim1->world->addParticle({ 3.8232564316 , 0.623837854536 , -2.81231187829 },0);
		sim2->world->addParticle({ 3.8232564316 , 0.623837854536 , -2.81231187829 },0);
		sim1->world->addParticle({ -3.69677943717 , -4.39807433257 , -2.89714251831 },0);
		sim2->world->addParticle({ -3.69677943717 , -4.39807433257 , -2.89714251831 },0);
		sim1->world->addParticle({ 2.08284715834 , 1.23535125907 , -0.939491972939 },0);
		sim2->world->addParticle({ 2.08284715834 , 1.23535125907 , -0.939491972939 },0);
		sim1->world->addParticle({ -1.20797796272 , -0.654588079093 , -0.558024570279 },0);
		sim2->world->addParticle({ -1.20797796272 , -0.654588079093 , -0.558024570279 },0);
		sim1->world->addParticle({ 3.23287437642 , 0.513286541268 , 2.9453553999 },0);
		sim2->world->addParticle({ 3.23287437642 , 0.513286541268 , 2.9453553999 },0);

		sim1->calculateInteractionForcesEnergiesWithLattice(4);
		sim2->calculateInteractionForcesEnergiesNaive();

		TS_ASSERT_DELTA(sim1->world->energy, sim2->world->energy, 1e-12);
		for (unsigned int i=0; i<64; i++) {
			TS_ASSERT_DELTA(
				sim1->world->activeParticles[i].cumulativeForce,
				sim2->world->activeParticles[i].cumulativeForce,
				1e-12);
		}

		delete sim1;
		delete sim2;
	}

	/* same as before but with LennardJones particles */
	void test_neighborLattice_ForcesEnergiesSameAsWithout_lj(void)
	{
		Simulation * sim1 = new Simulation();
		Simulation * sim2 = new Simulation();
		sim1->config->boxsize = 10.;
		sim2->config->boxsize = 10.;
		sim1->config->new_Type("A", .44, 1., 1.);
		sim2->config->new_Type("A", .44, 1., 1.);
		sim1->config->new_LennardJones("A<->A", {0,0}, 1.);
		sim2->config->new_LennardJones("A<->A", {0,0}, 1.);

		sim1->world->addParticle({ 4.76666031119 , 4.0160170777 , 4.0672375518 },0);
		sim2->world->addParticle({ 4.76666031119 , 4.0160170777 , 4.0672375518 },0);
		sim1->world->addParticle({ -3.98539898383 , 2.73922988841 , -1.88756863117 },0);
		sim2->world->addParticle({ -3.98539898383 , 2.73922988841 , -1.88756863117 },0);
		sim1->world->addParticle({ -4.64113537674 , -4.77370165644 , 1.49917576232 },0);
		sim2->world->addParticle({ -4.64113537674 , -4.77370165644 , 1.49917576232 },0);
		sim1->world->addParticle({ 4.35232059755 , -1.26274520297 , 3.44115107475 },0);
		sim2->world->addParticle({ 4.35232059755 , -1.26274520297 , 3.44115107475 },0);
		sim1->world->addParticle({ 3.68982485953 , -1.93423848154 , -4.14627025789 },0);
		sim2->world->addParticle({ 3.68982485953 , -1.93423848154 , -4.14627025789 },0);
		sim1->world->addParticle({ -4.79815051305 , -3.90603258185 , -0.097594026613 },0);
		sim2->world->addParticle({ -4.79815051305 , -3.90603258185 , -0.097594026613 },0);
		sim1->world->addParticle({ -4.88100872325 , -2.0740110649 , -3.6106471889 },0);
		sim2->world->addParticle({ -4.88100872325 , -2.0740110649 , -3.6106471889 },0);
		sim1->world->addParticle({ -4.42779240701 , -2.02237343148 , 1.51030607448 },0);
		sim2->world->addParticle({ -4.42779240701 , -2.02237343148 , 1.51030607448 },0);
		sim1->world->addParticle({ 3.47709914903 , -1.73181016779 , 4.52081071425 },0);
		sim2->world->addParticle({ 3.47709914903 , -1.73181016779 , 4.52081071425 },0);
		sim1->world->addParticle({ -0.389747230187 , -3.82980299407 , -2.49044443826 },0);
		sim2->world->addParticle({ -0.389747230187 , -3.82980299407 , -2.49044443826 },0);
		sim1->world->addParticle({ 3.59121711402 , -2.90430855842 , -0.14096088793 },0);
		sim2->world->addParticle({ 3.59121711402 , -2.90430855842 , -0.14096088793 },0);
		sim1->world->addParticle({ -3.55055171936 , 2.88152071785 , 4.45628407752 },0);
		sim2->world->addParticle({ -3.55055171936 , 2.88152071785 , 4.45628407752 },0);
		sim1->world->addParticle({ 3.24509853495 , 4.02072143176 , -0.373666874198 },0);
		sim2->world->addParticle({ 3.24509853495 , 4.02072143176 , -0.373666874198 },0);
		sim1->world->addParticle({ 3.88166132438 , 3.00271278363 , -1.17663972183 },0);
		sim2->world->addParticle({ 3.88166132438 , 3.00271278363 , -1.17663972183 },0);
		sim1->world->addParticle({ 2.54568685315 , 4.32619498531 , 3.06384782357 },0);
		sim2->world->addParticle({ 2.54568685315 , 4.32619498531 , 3.06384782357 },0);
		sim1->world->addParticle({ -2.88555772724 , 1.77393019522 , -1.26935471838 },0);
		sim2->world->addParticle({ -2.88555772724 , 1.77393019522 , -1.26935471838 },0);
		sim1->world->addParticle({ 3.77696938922 , 3.16848330108 , -4.37568115709 },0);
		sim2->world->addParticle({ 3.77696938922 , 3.16848330108 , -4.37568115709 },0);
		sim1->world->addParticle({ 3.43569668589 , 3.83653760004 , -3.23323736704 },0);
		sim2->world->addParticle({ 3.43569668589 , 3.83653760004 , -3.23323736704 },0);
		sim1->world->addParticle({ 3.40467837923 , -3.3114072018 , -1.88288324683 },0);
		sim2->world->addParticle({ 3.40467837923 , -3.3114072018 , -1.88288324683 },0);
		sim1->world->addParticle({ -3.64837816252 , -0.513877524659 , 3.44198156094 },0);
		sim2->world->addParticle({ -3.64837816252 , -0.513877524659 , 3.44198156094 },0);
		sim1->world->addParticle({ 1.61741126688 , -1.49326456399 , -4.0617649499 },0);
		sim2->world->addParticle({ 1.61741126688 , -1.49326456399 , -4.0617649499 },0);
		sim1->world->addParticle({ -3.32945494155 , -1.87248781176 , -2.68647512905 },0);
		sim2->world->addParticle({ -3.32945494155 , -1.87248781176 , -2.68647512905 },0);
		sim1->world->addParticle({ 1.689144414 , -3.51805286989 , -0.264814465122 },0);
		sim2->world->addParticle({ 1.689144414 , -3.51805286989 , -0.264814465122 },0);
		sim1->world->addParticle({ 4.91198059148 , -2.35888627297 , 2.27512953701 },0);
		sim2->world->addParticle({ 4.91198059148 , -2.35888627297 , 2.27512953701 },0);
		sim1->world->addParticle({ -0.813714232141 , 0.575708060963 , 4.68565057205 },0);
		sim2->world->addParticle({ -0.813714232141 , 0.575708060963 , 4.68565057205 },0);
		sim1->world->addParticle({ -0.282141268778 , 3.45753506633 , -2.60577987214 },0);
		sim2->world->addParticle({ -0.282141268778 , 3.45753506633 , -2.60577987214 },0);
		sim1->world->addParticle({ -3.98276394199 , 0.407239295958 , -2.7031966429 },0);
		sim2->world->addParticle({ -3.98276394199 , 0.407239295958 , -2.7031966429 },0);
		sim1->world->addParticle({ -2.52595699252 , -0.581450169204 , 3.6930847048 },0);
		sim2->world->addParticle({ -2.52595699252 , -0.581450169204 , 3.6930847048 },0);
		sim1->world->addParticle({ -1.11286690952 , 4.64944848343 , 4.11270017642 },0);
		sim2->world->addParticle({ -1.11286690952 , 4.64944848343 , 4.11270017642 },0);
		sim1->world->addParticle({ -3.18200236735 , 2.1999183697 , -4.44321771201 },0);
		sim2->world->addParticle({ -3.18200236735 , 2.1999183697 , -4.44321771201 },0);
		sim1->world->addParticle({ 4.50524340341 , 1.9473363275 , 1.5552971865 },0);
		sim2->world->addParticle({ 4.50524340341 , 1.9473363275 , 1.5552971865 },0);
		sim1->world->addParticle({ 0.259520230474 , 2.93920613198 , 1.04958655299 },0);
		sim2->world->addParticle({ 0.259520230474 , 2.93920613198 , 1.04958655299 },0);
		sim1->world->addParticle({ 2.68201041657 , -1.60438429561 , 3.85658334308 },0);
		sim2->world->addParticle({ 2.68201041657 , -1.60438429561 , 3.85658334308 },0);
		sim1->world->addParticle({ -4.53612481004 , -4.00778172216 , 1.04873965225 },0);
		sim2->world->addParticle({ -4.53612481004 , -4.00778172216 , 1.04873965225 },0);
		sim1->world->addParticle({ -1.34125683486 , 4.04621830104 , 2.61736270126 },0);
		sim2->world->addParticle({ -1.34125683486 , 4.04621830104 , 2.61736270126 },0);
		sim1->world->addParticle({ 3.64970978204 , -3.71023140909 , 1.14881718916 },0);
		sim2->world->addParticle({ 3.64970978204 , -3.71023140909 , 1.14881718916 },0);
		sim1->world->addParticle({ 0.951420026083 , -3.60219562428 , 4.56721030973 },0);
		sim2->world->addParticle({ 0.951420026083 , -3.60219562428 , 4.56721030973 },0);
		sim1->world->addParticle({ 0.566532770478 , -3.59371976143 , -0.399224116342 },0);
		sim2->world->addParticle({ 0.566532770478 , -3.59371976143 , -0.399224116342 },0);
		sim1->world->addParticle({ 0.643601927649 , -1.63725014476 , -2.55780192644 },0);
		sim2->world->addParticle({ 0.643601927649 , -1.63725014476 , -2.55780192644 },0);
		sim1->world->addParticle({ -1.16929995696 , -1.50883600191 , -4.60613066258 },0);
		sim2->world->addParticle({ -1.16929995696 , -1.50883600191 , -4.60613066258 },0);
		sim1->world->addParticle({ 0.115257329108 , 0.0971088102476 , 3.82966022532 },0);
		sim2->world->addParticle({ 0.115257329108 , 0.0971088102476 , 3.82966022532 },0);
		sim1->world->addParticle({ -0.819817935497 , 1.24197319149 , -3.45945437102 },0);
		sim2->world->addParticle({ -0.819817935497 , 1.24197319149 , -3.45945437102 },0);
		sim1->world->addParticle({ 0.514532520888 , -0.72326440747 , -1.09328218488 },0);
		sim2->world->addParticle({ 0.514532520888 , -0.72326440747 , -1.09328218488 },0);
		sim1->world->addParticle({ -1.64474217114 , 2.59855339803 , 3.20549028185 },0);
		sim2->world->addParticle({ -1.64474217114 , 2.59855339803 , 3.20549028185 },0);
		sim1->world->addParticle({ -4.52534872793 , 2.83906426899 , -2.78964619962 },0);
		sim2->world->addParticle({ -4.52534872793 , 2.83906426899 , -2.78964619962 },0);
		sim1->world->addParticle({ 1.26422676362 , -2.50418742356 , 0.741772054412 },0);
		sim2->world->addParticle({ 1.26422676362 , -2.50418742356 , 0.741772054412 },0);
		sim1->world->addParticle({ 0.0896373137543 , 2.51272514604 , -2.27164583661 },0);
		sim2->world->addParticle({ 0.0896373137543 , 2.51272514604 , -2.27164583661 },0);
		sim1->world->addParticle({ 2.31052356364 , 4.03974594927 , 0.180380518953 },0);
		sim2->world->addParticle({ 2.31052356364 , 4.03974594927 , 0.180380518953 },0);
		sim1->world->addParticle({ -4.78240792721 , -4.47028600724 , 2.44929179718 },0);
		sim2->world->addParticle({ -4.78240792721 , -4.47028600724 , 2.44929179718 },0);
		sim1->world->addParticle({ 4.39062304743 , 2.12872471111 , -3.74572811228 },0);
		sim2->world->addParticle({ 4.39062304743 , 2.12872471111 , -3.74572811228 },0);
		sim1->world->addParticle({ 4.2077143449 , 4.38691778031 , 0.389129309429 },0);
		sim2->world->addParticle({ 4.2077143449 , 4.38691778031 , 0.389129309429 },0);
		sim1->world->addParticle({ 1.04483848843 , -4.0409645901 , 1.58718447085 },0);
		sim2->world->addParticle({ 1.04483848843 , -4.0409645901 , 1.58718447085 },0);
		sim1->world->addParticle({ 1.76733178783 , 0.703987051801 , 3.87663701145 },0);
		sim2->world->addParticle({ 1.76733178783 , 0.703987051801 , 3.87663701145 },0);
		sim1->world->addParticle({ 4.59309977994 , -1.79952661068 , 1.22663557472 },0);
		sim2->world->addParticle({ 4.59309977994 , -1.79952661068 , 1.22663557472 },0);
		sim1->world->addParticle({ 2.78734997525 , -2.02091109424 , 0.545112115803 },0);
		sim2->world->addParticle({ 2.78734997525 , -2.02091109424 , 0.545112115803 },0);
		sim1->world->addParticle({ 2.05717043923 , -3.31714421527 , 2.6164577159 },0);
		sim2->world->addParticle({ 2.05717043923 , -3.31714421527 , 2.6164577159 },0);
		sim1->world->addParticle({ 2.53709775181 , -1.98531757353 , 4.75641140986 },0);
		sim2->world->addParticle({ 2.53709775181 , -1.98531757353 , 4.75641140986 },0);
		sim1->world->addParticle({ 2.72627079184 , -2.20004768919 , -3.11068702192 },0);
		sim2->world->addParticle({ 2.72627079184 , -2.20004768919 , -3.11068702192 },0);
		sim1->world->addParticle({ 3.79352596159 , -4.95754411006 , -1.47110394297 },0);
		sim2->world->addParticle({ 3.79352596159 , -4.95754411006 , -1.47110394297 },0);
		sim1->world->addParticle({ 0.239808202985 , 2.0312977463 , 0.0521382879886 },0);
		sim2->world->addParticle({ 0.239808202985 , 2.0312977463 , 0.0521382879886 },0);
		sim1->world->addParticle({ 1.12660263915 , 2.89963130119 , -2.77052361871 },0);
		sim2->world->addParticle({ 1.12660263915 , 2.89963130119 , -2.77052361871 },0);
		sim1->world->addParticle({ -4.53936736303 , 4.51646850838 , 0.74203346798 },0);
		sim2->world->addParticle({ -4.53936736303 , 4.51646850838 , 0.74203346798 },0);
		sim1->world->addParticle({ 3.34606284337 , 4.95847490794 , -0.649037847239 },0);
		sim2->world->addParticle({ 3.34606284337 , 4.95847490794 , -0.649037847239 },0);
		sim1->world->addParticle({ 1.46424626004 , 2.87204978011 , 2.29630572042 },0);
		sim2->world->addParticle({ 1.46424626004 , 2.87204978011 , 2.29630572042 },0);
		
		sim1->calculateInteractionForcesEnergiesWithLattice(4);
		sim2->calculateInteractionForcesEnergiesNaive();

		TS_ASSERT_DELTA(sim1->world->energy, sim2->world->energy, 1e-12);
		for (unsigned int i=0; i<64; i++) {
			TS_ASSERT_DELTA(
				sim1->world->activeParticles[i].cumulativeForce,
				sim2->world->activeParticles[i].cumulativeForce,
				1e-12);
		}
		
		delete sim1;
		delete sim2;
	}

};
