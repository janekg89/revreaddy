/* Fusion3.h  A+B <-> C
 * A reaction that is bimolecular in forward direction
 * and unimolecular in backward direction. */

#ifndef __FUSION3_H_INCLUDED__
#define __FUSION3_H_INCLUDED__
#include "Reaction.h"
#include <vector>
#include <memory>
#include <string>
#include "ParticleInteraction.h"

class Fusion3 : public Reaction
{
	public:
	Fusion3(
		std::string inName,
		std::vector<unsigned> inForwardTypes,
		std::vector<unsigned> inBackwardTypes,
		double inForwardRate,
		double inBackwardRate,
		Random * inRandom);

	/* pointers to the interactions between particles A and B */
	std::vector< std::shared_ptr<ParticleInteraction> > interactions;
	/* inverse partition function value for drawing the distance
	 * of A and B particle after backward reaction, must be 
	 * correctly set according to the interaction between A and B */
	double inversePartition;
	/* The maximum of the probability distribution is needed for
	 * drawing from that distribution. Must be set correctly 
	 * from the outside */
	double maxDistr;
	/* Sum of the radii of particles A and B */
	double radiiSum;
	/* Sum of the reactionRadii of particles A and B */
	double reactionRadiiSum;
	/* The mean value of the distribution is returned if no
	 * random value could be successfully drawn */
	double meanDistr;
	/* 1 / k_BT */
	double inverseTemperature;
	/* if true, an interaction has been assigned (not 
	 * guaranteeing if it's the only one possible) and 
	 * the other parameters have been set */
	bool isConfigured;
	/* The particles radii to place new particles at center of mass.
	 * The weights correspond to the particles masses */
	double radiusA;
	double radiusB;
	double weightA;
	double weightB;
	/* configure() sets the parameters that depend on
	 * energy functions of particle interactions */
	void configure(
		std::vector< std::shared_ptr<ParticleInteraction> > inInteractions,	
		double inInversePartition,
		double inMaxDistr,
		double inRadiiSum,
		double inReactionRadiiSum,
		double inMeanDistr,
		double inInverseTemperature,
		double inRadiusA,
		double inRadiusB);
	double performForward(
		std::vector<unsigned long> particleIndices,
		World * world,
		double timestep);
	double performBackward(
		std::vector<unsigned long> particleIndices,
		World * world,
		double timestep);

	private:
	/* f(x) = Z^-1 * x * exp[ -beta * potential( x * radiiSum ) ] */
	double distribution(double x);
	/* draw a uniform number from the distribution above in the range [0,1] */
	double uniformFromDistribution();
};

#endif // __FUSION3_H_INCLUDED__