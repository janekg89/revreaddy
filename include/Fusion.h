/* Fusion.h  A+B <-> C
 * A reaction that is bimolecular in forward direction
 * and unimolecular in backward direction. */

#ifndef __FUSION_H_INCLUDED__
#define __FUSION_H_INCLUDED__
#include "Reaction.h"
#include <vector>
#include <string>
#include "ParticleInteraction.h"

class Fusion : public Reaction
{
	public:
	Fusion(
		std::string inName,
		std::vector<unsigned> inForwardTypes,
		std::vector<unsigned> inBackwardTypes,
		double inForwardRate,
		double inBackwardRate,
		Random * inRandom);
	~Fusion();

	/* pointers to the interactions between particles A and B */
	std::vector<ParticleInteraction*> interactions;
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
	/* The mean value of the distribution is returned if no
	 * random value could be successfully drawn */
	double meanDistr;
	/* 1 / k_BT */
	double inverseTemperature;
	/* if true, an interaction has been assigned (not 
	 * guaranteeing if it's the only one possible) and 
	 * the other parameters have been set */
	bool isConfigured;
	/* configure() sets the parameters that depend on
	 * energy functions of particle interactions */
	void configure(
		std::vector<ParticleInteraction*> inInteractions,	
		double inInversePartition,
		double inMaxDistr,
		double inRadiiSum,
		double inMeanDistr,
		double inInverseTemperature);
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

#endif // __FUSION_H_INCLUDED__