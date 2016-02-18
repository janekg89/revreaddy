/* Random.cpp */
#include "Random.h"

Random::Random(std::string type) {
	gsl_rng_env_setup();

	if (type == "mt19937")	{ this->randomGeneratorType = gsl_rng_mt19937; }
	else if (type == "taus") { this->randomGeneratorType = gsl_rng_taus; }
	else if (type == "ranlxs0") { this->randomGeneratorType = gsl_rng_ranlxs0; }
	else { this->randomGeneratorType = gsl_rng_ranlxs0; }

	this->randomGeneratorHandle = gsl_rng_alloc(this->randomGeneratorType);
	this->toSeed();
	gsl_rng_set(this->randomGeneratorHandle, this->seed);
	LOG_INFO(
		"Random number generator of type '" 
		<< gsl_rng_name(randomGeneratorHandle) << "' has been intialized and seeded."
	)
}

Random::~Random() {
	gsl_rng_free(this->randomGeneratorHandle);
	LOG_INFO("Random number generator has been released.")
}

void Random::toSeed() {
	FILE *devurandom;
	devurandom = fopen("/dev/urandom","r");
	fread(&this->seed, sizeof(this->seed), 1, devurandom);
	fclose(devurandom);
	LOG_INFO("Got seed " << this->seed << " from /dev/urandom.")
}

double Random::normal() {
	return gsl_ran_gaussian_ziggurat(this->randomGeneratorHandle, 1.0);
}

std::vector<double> Random::normal3D() {
	std::vector<double> randomArray = {0.,0.,0.};
	randomArray[0] = gsl_ran_gaussian_ziggurat(this->randomGeneratorHandle,1.0);
	randomArray[1] = gsl_ran_gaussian_ziggurat(this->randomGeneratorHandle,1.0);
	randomArray[2] = gsl_ran_gaussian_ziggurat(this->randomGeneratorHandle,1.0);
	return randomArray;
}

double Random::uniform() {
	return gsl_rng_uniform(this->randomGeneratorHandle);
}
