/* Random.h
 * author: Christoph Froehner
 *
 * This object will construct a random number generator,
 * from which normal distributed and uniformly distributed
 * random numbers/arrays are drawn.
 * At construction the random generator is seeded by
 * /dev/urandom.
 *
 */

#ifndef __RANDOM_H_INCLUDED__
#define __RANDOM_H_INCLUDED__
#include <gsl/gsl_rng.h> 	// GNU random number generators
#include <gsl/gsl_randist.h>// GNU random number distributions
#include <string>
#include <array>

class Random
{
	private:
		const gsl_rng_type * randomGeneratorType;
		gsl_rng * randomGeneratorHandle;
		unsigned long int seed;
	public:
		// type can be "mt19937", "taus" or "ranlxs0"
		Random(std::string type);
		~Random();
		void		toSeed();
		double		normal();
		std::array<double, 3> normal3D();
		double 		uniform();
};
#endif
