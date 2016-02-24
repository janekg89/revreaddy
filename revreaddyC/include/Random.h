/* Random.h
 * This object will construct a random number generator,
 * from which normal distributed and uniformly distributed
 * random numbers/arrays are drawn.
 * At construction the random generator is seeded by
 * /dev/urandom on unix systems. */

#ifndef __RANDOM_H_INCLUDED__
#define __RANDOM_H_INCLUDED__
#include <string>
#include <vector>
#include <iostream>
#include <gsl/gsl_rng.h> // GNU random number generators
#include <gsl/gsl_randist.h> // GNU random number distributions
#include "logging.h"

class Random {
private:
	const gsl_rng_type * randomGeneratorType;
	gsl_rng * randomGeneratorHandle;
	unsigned long int seed;
public:
	// type can be "mt19937", "taus" or "ranlxs0"
	Random(std::string type = "ranlxs0");
	~Random();
	std::string getType();
	unsigned long int getNewSeed();
	double normal();
	std::vector<double> normal3D();
	double uniform();
};

#endif // __RANDOM_H_INCLUDED__