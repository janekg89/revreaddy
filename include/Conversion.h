/* Conversion.h
 * An implementation of Reaction, which is unimolecular
 * in both directions. The affected particle is destroyed
 * and at the same position the new particle is created. */

#ifndef __CONVERSION_H_INCLUDED__
#define __CONVERSION_H_INCLUDED__
#include "Reaction.h"
#include <vector>
#include <string>

class Conversion : public Reaction
{
	public:
		Conversion(
			std::string inName,
			std::vector<unsigned int> inForwardTypes,
			std::vector<unsigned int> inBackwardTypes,
			double inForwardRate,
			double inBackwardRate,
			Random * inRandom);
		~Conversion();

		double performForward(
			std::vector<unsigned long int> particleIndices,
			World * world,
			double timestep);
		double performBackward(
			std::vector<unsigned long int> particleIndices,
			World * world,
			double timestep);
};

#endif //__CONVERSION_H_INCLUDED__