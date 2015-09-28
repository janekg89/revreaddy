/* Config.cpp */

#include "Config.h"

Config::Config(World * inWorld, Random * inRandom)
{
	this->world                 = inWorld;
	this->random                = inRandom;
	this->timestep              = 0.001;
	this->temperature           = 1.;
	this->kBoltzmann            = 1.;
	this->isPeriodic            = true;
	this->boxsize               = 10.;
	this->isReversibleDynamics  = true;
	this->isReversibleReactions = true;
	this->useNeighborList       = true;
	this->reactionPropagation   = 0;
}

Config::~Config()
{
	this->deleteAllObservables();
	this->deleteAllGeometries();
	this->deleteAllForces();
	this->deleteAllReactions();
}

void Config::new_Type(
	std::string name,
	double radius,
	double diffusionConstant,
	double reactionRadius)
{
	if ( radius < 0. ) {
		std::cout << "Error: The particle radius must be non-negative!\n"
		          << "New type is not created" << std::endl;
		return;
	}
	if ( diffusionConstant < 0. ) {
		std::cout << "Error: The diffusionConstant must be non-negative!\n"
		          << "New Type is not created" << std::endl;
		return;
	}
	if ( reactionRadius < 0. ) {
		std::cout << "Error: The reactionRadius must be non-negative!\n"
		          << "New Type is not created" << std::endl;
		return;
	}
	ParticleType pType(
		name,
		radius,
		diffusionConstant,
		reactionRadius);
	this->typeDict.push_back(pType);
}

unsigned int Config::getNumberOfTypes() {
	return this->typeDict.size();
}

std::string Config::getDictName(unsigned int i) {
	return this->typeDict[i].name;
}

double Config::getDictRadius(unsigned int i) {
	return this->typeDict[i].radius;
}

double Config::getDictDiffusionConstant(unsigned int i) {
	return this->typeDict[i].diffusionConstant;
}

double Config::getDictReactionRadius(unsigned int i) {
	return this->typeDict[i].reactionRadius;
}

unsigned int Config::getParticleNumber() {
	return world->activeParticles.size();
}


void Config::writeAllObservablesToFile()
{
	for (auto* obs : this->observables) {
		obs->writeBufferToFile();
	}
}

void Config::writeLastObservableToFile()
{
	if (this->observables.size() > 0) {
		this->observables.back()->writeBufferToFile();
	}
	else {
		std::cout << "Error: There are no observables to write" << std::endl;
	}
}

std::string Config::showObservables()
{
	std::string content = "Observables: ";
	if (this->observables.size() > 0) {
		for (auto* obs : this->observables) {
			content += std::string (typeid(*obs).name()) + " ";
		}
	}
	else {content += "empty";}
	content += "\n";
	return content;
}

void Config::deleteAllObservables()
{
	/* first delete the obs, since they we're allocated with 'new'
	 * then erase the pointers from the vector */
	for (auto* obs : this->observables) {
		delete obs;
	}
	this->observables.clear();
}

void Config::deleteLastObservable()
{
	/* first delete the observable, then its pointer in the vector */
	if (this->observables.size() > 0) {
		delete this->observables.back();
		this->observables.pop_back();
	}
}

void Config::new_Trajectory(unsigned long int recPeriod, std::string filename)
{
	Trajectory * obs = new Trajectory(filename);
	obs->recPeriod = recPeriod;
	this->observables.push_back(obs);
}

void Config::new_RadialDistribution(
	unsigned long int recPeriod,
	std::string filename,
	std::vector<double> ranges,
	std::vector< std::vector<unsigned int> > considered)
{
	RadialDistribution * rad = new RadialDistribution(
		ranges,
		this->isPeriodic,
		this->boxsize,
		considered,
		filename);
	rad->recPeriod = recPeriod;
	this->observables.push_back(rad);
}

void Config::new_MeanSquaredDisplacement(
	unsigned long int recPeriod,
	std::string filename,
	unsigned int particleTypeId)
{
	MeanSquaredDisplacement * msd = new MeanSquaredDisplacement(
		world->activeParticles,
		particleTypeId,
		world->cumulativeRuntime,
		this->boxsize,
		filename);
	msd->recPeriod = recPeriod;
	this->observables.push_back(msd);
}

void Config::new_ProbabilityDensity(
	unsigned long recPeriod,
	std::string filename,
	unsigned particleTypeId,
	std::vector<double> range,
	unsigned coord)
{
	ProbabilityDensity * prob = new ProbabilityDensity(
		recPeriod,
		0,
		filename,
		particleTypeId,
		range,
		coord,
		this->world);
	this->observables.push_back(prob);
}

void Config::new_Energy(unsigned long int recPeriod, std::string filename)
{
	Energy * ener = new Energy(
		recPeriod,
		0,
		filename);
	this->observables.push_back(ener);
}

void Config::new_Acceptance(unsigned long int recPeriod, std::string filename)
{
	Acceptance * acc = new Acceptance(
		recPeriod,
		0,
		filename);
	this->observables.push_back(acc);
}

void Config::new_ParticleNumbers(
	unsigned long recPeriod,
	std::string filename,
	unsigned particleTypeId)
{
	ParticleNumbers * par = new ParticleNumbers(
		recPeriod,
		0,
		filename,
		particleTypeId);
	this->observables.push_back(par);
}

void Config::deleteAllGeometries()
{
	/* first delete the geometries, since they we're allocated with 'new'
	 * then erase the pointers from the vector */
	for (auto* geo : this->geometries) {
		delete geo;
	}
	this->geometries.clear();
}

void Config::new_Wall(
	std::vector<double> normal,
	std::vector<double> point,
	double strength,
	std::vector<unsigned int> particleTypeIds)
{
	Wall * wall = new Wall(
		normal,
		point,
		strength,
		particleTypeIds);
	this->geometries.push_back(wall);
}

void Config::new_DoubleWellZ(
	double distanceMinima,
	double strength,
	std::vector<unsigned int> particleTypeIds)
{
	DoubleWellZ * well = new DoubleWellZ(
		distanceMinima,
		strength,
		particleTypeIds);
	this->geometries.push_back(well);
}

void Config::deleteAllForces()
{
	/* first delete the forces, since they we're allocated with 'new'
	 * then erase the pointers from the vector */
	for (auto* f : this->possibleInteractions) {
		delete f;
	}
	this->possibleInteractions.clear();
}

void Config::new_SoftRepulsion(
	std::string name,
	std::vector<unsigned int> affectedTuple,
	double repulsionStrength)
{
	if (affectedTuple.size() != 2) {
		std::cout << "Error: The given tuple must be of length 2" << std::endl;
		return;
	}
	if ( (affectedTuple[0] > ( this->typeDict.size() - 1) ) 
	  || (affectedTuple[1] > ( this->typeDict.size() - 1) ) ) {
		std::cout << "Error: The given particle type(s) do not exist.\n"
		          << "Make sure to add them first" << std::endl;
		return;
	}
	if ( repulsionStrength <= 0. ) {
		std::cout << "Error: The repulsion strength must be larger than zero"
		          << std::endl;
		return;
	}
	SoftRepulsion * soft = new SoftRepulsion(
		name,
		affectedTuple,
		repulsionStrength);
	// set cutoff correctly
	soft->cutoff = this->typeDict[affectedTuple[0]].radius 
	             + this->typeDict[affectedTuple[1]].radius;
	this->possibleInteractions.push_back(soft);
	std::cout << "Info: SoftRepulsion interaction added to possibleInteractions"
	          << std::endl;
}

void Config::new_LennardJones(
	std::string name,
	std::vector<unsigned int> affectedTuple,
	double epsilon)
{
	if (affectedTuple.size() != 2) {
		std::cout << "Error: The given tuple must be of length 2" << std::endl;
		return;
	}
	if ( (affectedTuple[0] > ( this->typeDict.size() - 1) ) 
	  || (affectedTuple[1] > ( this->typeDict.size() - 1) ) ) {
		std::cout << "Error: The given particle type(s) do not exist.\n"
		          << "Make sure to add them first" << std::endl;
		return;
	}
	if ( epsilon <= 0. ) {
		std::cout << "Error: The given epsilon must be larger than zero"
		          << std::endl;
		return;
	}
	LennardJones * lj = new LennardJones(
		name,
		affectedTuple,
		epsilon);
	// set cutoff correctly
	lj->cutoff = 2.5 * ( this->typeDict[affectedTuple[0]].radius
	                   + this->typeDict[affectedTuple[1]].radius );
	this->possibleInteractions.push_back(lj);
	std::cout << "Info: LennardJones interaction added to possibleInteractions"
	          << std::endl;
}

unsigned int Config::getNumberForces()
{
	return this->possibleInteractions.size();
}

std::string Config::getForceName(unsigned int i)
{
	return this->possibleInteractions[i]->name;
}

std::string Config::getForceType(unsigned int i)
{
	return this->possibleInteractions[i]->type;
}

std::vector<unsigned int> Config::getForceAffectedTuple(unsigned int i)
{
	return this->possibleInteractions[i]->affectedTuple;
}

std::vector<double> Config::getForceParameters(unsigned int i)
{
	return this->possibleInteractions[i]->parameters;
}

void Config::deleteAllReactions()
{
	/* first delete the reactions, since they we're allocated with 'new'
	 * then erase the pointers from the vector */
	for (auto* reaction : this->possibleReactions) {
		delete reaction;
	}
	this->possibleReactions.clear();	
}

void Config::new_Conversion(
	std::string name,
	unsigned forwardType,
	unsigned backwardType,
	double forwardRate,
	double backwardRate)
{
	std::vector<unsigned> forwardTypes;
	forwardTypes.push_back(forwardType);
	std::vector<unsigned> backwardTypes;
	backwardTypes.push_back(backwardType);
	Conversion * conv = new Conversion(
		name,
		forwardTypes,
		backwardTypes,
		forwardRate,
		backwardRate,
		this->random);
	this->possibleReactions.push_back(conv);
}