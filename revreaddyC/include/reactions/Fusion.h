/* Fusion.h  A+B <-> C
 * A reaction that is bimolecular in forward direction
 * and unimolecular in backward direction. */

#ifndef __FUSION_H_INCLUDED__
#define __FUSION_H_INCLUDED__
#include <vector>
#include <memory>
#include <string>
#include "Reaction.h"
#include "World.h"
#include "Random.h"
#include "ParticleInteraction.h"
#include "utils.h"
#include "logging.h"

class Fusion : public Reaction {
public:
	Fusion(
		std::string inName,
		std::vector<unsigned> inForwardTypes,
		std::vector<unsigned> inBackwardTypes,
		double inForwardRate,
		double inBackwardRate,
		double inReactionDistance);
	~Fusion();

	/* pointers to the interactions between particles A and B */
	std::vector< std::shared_ptr<ParticleInteraction> > interactions;
	/* pointer to utilities needed for minimum distance */
	Utils * utils;
	/* inverse partition function value for drawing the distance
	 * of A and B particle after backward reaction, must be 
	 * correctly set according to the interaction between A and B */
	double inversePartition;
	/* The maximum of the probability distribution is needed for
	 * drawing from that distribution. Must be set correctly 
	 * from the outside */
	double maxDistr;
	/* The mean value of the distribution is returned if no
	 * random value could be successfully drawn */
	double meanDistr;
	/* 1 / k_BT */
	double inverseTemperature;
	/* The particles radii to place new particles at center of mass.
	 * The weights correspond to the particles masses */
	double radiusA;
	double radiusB;
	double weightA;
	double weightB;
    double radiiSum;
	bool isPeriodic;
	double boxsize;
	/* configure() sets the parameters that depend on
	 * energy functions of particle interactions */
	void configure(
		std::vector< std::shared_ptr<ParticleInteraction> > inInteractions,	
		double inInversePartition,
		double inMaxDistr,
		double inMeanDistr,
		double inInverseTemperature,
		double inRadiusA,
		double inRadiusB,
		bool inIsPeriodic,
		double inBoxsize);
	double performForward(
		std::vector<unsigned long> particleIndices,
		double timestep,
		World * world,
		Random * random);
	double performBackward(
		std::vector<unsigned long> particleIndices,
		double timestep,
		World * world,
		Random * random);

private:
	/* f(x) = Z^-1 * exp[ -beta * potential( x ) ] */
	double distribution(double x);
	/* draw a number from the distribution above in [0,reactionRadiiSum] */
	double randomFromDistribution(Random * random);
};

#endif // __FUSION_H_INCLUDED__