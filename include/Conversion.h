/* Conversion.h
 * An implementation of Reaction, which is unimolecular
 * in both directions. */

#ifndef __CONVERSION_H_INCLUDED__
#define __CONVERSION_H_INCLUDED__
#include "Reaction.h"

class Conversion : public Reaction
{
	public:
		Conversion(
			std::string inName,
			std::vector<unsigned int> inForwardTypes,
			std::vector<unsigned int> inBackwardTypes,
			double inForwardRate,
			double inBackwardRate);
		
		double performForward(
			std::vector<unsigned long long> particleIndices,
			Simulation* simulation);
		double performBackward(
			std::vector<unsigned long long> particleIndices,
			Simulation* simulation);
};

#endif //__CONVERSION_H_INCLUDED__
