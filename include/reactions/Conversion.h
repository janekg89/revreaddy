/* Conversion.h  A <-> B
 * An implementation of Reaction, which is unimolecular
 * in both directions. The affected particle is destroyed
 * and at the same position the new particle is created. */

#ifndef __CONVERSION_H_INCLUDED__
#define __CONVERSION_H_INCLUDED__
#include <vector>
#include <string>
#include "Reaction.h"
#include "World.h"
#include "Random.h"

class Conversion : public Reaction
{
public:
	Conversion(
		std::string inName,
		std::vector<unsigned int> inForwardTypes,
		std::vector<unsigned int> inBackwardTypes,
		double inForwardRate,
		double inBackwardRate);
	~Conversion();

	double performForward(
		std::vector<unsigned long int> particleIndices,
		double timestep,
		World * world,
		Random * random);
	double performBackward(
		std::vector<unsigned long int> particleIndices,
		double timestep,
		World * world,
		Random * random);
};

#endif //__CONVERSION_H_INCLUDED__