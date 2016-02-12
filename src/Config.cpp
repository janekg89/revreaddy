/* Config.cpp */

#include "Config.h"
#define print(x) std::cout << x << std::endl;

template<typename T, typename ...Args>
std::unique_ptr<T> make_unique( Args&& ...args )
{
	return std::unique_ptr<T>( new T( std::forward<Args>(args)... ) );
}

Config::Config()
{
	this->timestep              = 0.001;
	this->kT                    = 1.;
	this->isPeriodic            = true;
	this->boxsize               = 10.;
}

Config::~Config()
{
	this->deleteAllGeometries();
	this->deleteAllReactions();
	this->deleteAllInteractions();
}

void Config::new_Type(
	const std::string name,
	const double radius,
	const double diffusionConstant)
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
	ParticleType pType(
		name,
		radius,
		diffusionConstant);
	this->particleTypes.push_back(pType);
}

unsigned int Config::getNumberOfTypes() {
	return this->particleTypes.size();
}

std::string Config::getDictName(unsigned int i) {
	return this->particleTypes[i].name;
}

double Config::getDictRadius(unsigned int i) {
	return this->particleTypes[i].radius;
}

double Config::getDictDiffusionConstant(unsigned int i) {
	return this->particleTypes[i].diffusionConstant;
}

void Config::deleteAllGeometries()
{
	/* Erase all the unique pointers, the geometries are thus deleted as well */
	this->geometries.clear();
}

void Config::new_Wall(
	std::vector<double> normal,
	std::vector<double> point,
	double strength,
	std::vector<unsigned int> particleTypeIds)
{
	std::unique_ptr<Wall> wall = make_unique<Wall>(
		normal,
		point,
		strength,
		particleTypeIds);
	this->geometries.push_back( std::move(wall) );
}

void Config::new_DoubleWellZ(
	double distanceMinima,
	double strength,
	std::vector<unsigned int> particleTypeIds)
{
	std::unique_ptr<DoubleWellZ> well = make_unique<DoubleWellZ>(
		distanceMinima,
		strength,
		particleTypeIds);
	this->geometries.push_back( std::move(well) );
}

void Config::deleteAllInteractions()
{
	/* Erase all the shared pointers, if some interactions are still 
	 * referenced by Reactions, the interaction is technically valid,
	 * but it will lead to unwanted behavior. Therefore delete all
	 * reactions here as well. */
	this->reactions.clear();
	this->interactions.clear();
}

void Config::new_SoftRepulsion(
	std::string name,
	std::vector<unsigned> affectedTuple,
	double repulsionStrength)
{
	if (affectedTuple.size() != 2) {
		std::cout << "Error: The given tuple must be of length 2" << std::endl;
		return;
	}
	if ( (affectedTuple[0] > ( this->particleTypes.size() - 1) ) 
	  || (affectedTuple[1] > ( this->particleTypes.size() - 1) ) ) {
		std::cout << "Error: The given particle type(s) do not exist.\n"
		          << "Make sure to add them first" << std::endl;
		return;
	}
	if ( repulsionStrength <= 0. ) {
		std::cout << "Error: The repulsion strength must be larger than zero"
		          << std::endl;
		return;
	}
	std::shared_ptr<SoftRepulsion> soft = std::make_shared<SoftRepulsion>(
		name,
		affectedTuple,
		repulsionStrength);
	// set cutoff correctly
	soft->cutoff = this->particleTypes[affectedTuple[0]].radius 
	             + this->particleTypes[affectedTuple[1]].radius;
	this->interactions.push_back( std::move(soft) );
	std::cout << "Info: SoftRepulsion interaction added to interactions"
	          << std::endl;
}

void Config::new_LennardJones(
	std::string name,
	std::vector<unsigned> affectedTuple,
	double epsilon)
{
	if (affectedTuple.size() != 2) {
		std::cout << "Error: The given tuple must be of length 2" << std::endl;
		return;
	}
	if ( (affectedTuple[0] > ( this->particleTypes.size() - 1) ) 
	  || (affectedTuple[1] > ( this->particleTypes.size() - 1) ) ) {
		std::cout << "Error: The given particle type(s) do not exist.\n"
		          << "Make sure to add them first" << std::endl;
		return;
	}
	if ( epsilon <= 0. ) {
		std::cout << "Error: The given epsilon must be larger than zero"
		          << std::endl;
		return;
	}
	std::shared_ptr<LennardJones> lj = std::make_shared<LennardJones>(
		name,
		affectedTuple,
		epsilon);
	// set cutoff correctly
	lj->cutoff = 2.5 * ( this->particleTypes[affectedTuple[0]].radius
	                   + this->particleTypes[affectedTuple[1]].radius );
	this->interactions.push_back( std::move(lj) );
	std::cout << "Info: LennardJones interaction added to interactions"
	          << std::endl;
}

unsigned int Config::getNumberInteractions()
{
	return this->interactions.size();
}

std::string Config::getInteractionName(unsigned i)
{
	return this->interactions[i]->name;
}

std::string Config::getInteractionType(unsigned i)
{
	return this->interactions[i]->type;
}

std::vector<unsigned int> Config::getInteractionAffectedTuple(unsigned i)
{
	return this->interactions[i]->affectedTuple;
}

std::vector<double> Config::getInteractionParameters(unsigned i)
{
	return this->interactions[i]->parameters;
}

double Config::getInteractionCutoff(unsigned i)
{
	return this->interactions[i]->cutoff;
}

void Config::deleteAllReactions()
{
	/* Erase all the unique pointers from reactions. Thus all
	 * reactions will be destroyed accordingly */
	this->reactions.clear();	
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
	std::unique_ptr<Conversion> conv = make_unique<Conversion>(
		name,
		forwardTypes,
		backwardTypes,
		forwardRate,
		backwardRate);
	this->reactions.push_back( std::move(conv) );
}

void Config::new_Fusion3(
	std::string name,
	unsigned forwardTypeA,
	unsigned forwardTypeB,
	unsigned backwardTypeC,
	double forwardRate,
	double backwardRate,
	double reactionDistance)
{
	std::vector<unsigned> forwardTypes;
	forwardTypes.push_back(forwardTypeA);
	forwardTypes.push_back(forwardTypeB);
	std::vector<unsigned> backwardTypes;
	backwardTypes.push_back(backwardTypeC);
	std::unique_ptr<Fusion3> fus = make_unique<Fusion3>(
		name,
		forwardTypes,
		backwardTypes,
		forwardRate,
		backwardRate,
		reactionDistance);
	this->reactions.push_back( std::move(fus) );
}

/* That reactionIndex is of type Fusion3 has to be ensured by the caller */
void Config::configure_Fusion3(
	unsigned reactionIndex,
	std::vector<unsigned> interactionsIndices,
	double inversePartition,
	double maxDistr,
	double radiiSum,
	double reactionRadiiSum,
	double meanDistr,
	double inverseTemperature,
	double radiusA,
	double radiusB)
{
	std::vector< std::shared_ptr<ParticleInteraction> > consideredInteractions;
	for (unsigned i=0; i<interactionsIndices.size(); i++) {
		/* push_back constructs a new shared ptr to be
		 * forwarded to the Fusion3. */
		consideredInteractions.push_back(
			this->interactions[ interactionsIndices[i] ]
		);
	}
	// TODO this still seems inappropriate
	Fusion3 * fus = dynamic_cast<Fusion3*>( this->reactions[reactionIndex].get() );
	fus->configure(
		consideredInteractions,
		inversePartition,
		maxDistr,
		radiiSum,
		reactionRadiiSum,
		meanDistr,
		inverseTemperature,
		radiusA,
		radiusB,
		this->isPeriodic,
		this->boxsize);	
	/* Don't delete fus, since it is uniquely owned by Config. 
	 * The reference was only borrowed to cast it to Fusion3 
	 * and configure it properly */
}

unsigned Config::getNumberReactions()
{
	return this->reactions.size();
}

std::string Config::getReactionName(unsigned i)
{
	return this->reactions[i]->name;
}

std::string Config::getReactionType(unsigned i)
{
	return this->reactions[i]->type;
}

std::vector<unsigned> Config::getReactionForwardTypes(unsigned i)
{
	return this->reactions[i]->forwardTypes;
}

std::vector<unsigned> Config::getReactionBackwardTypes(unsigned i)
{
	return this->reactions[i]->backwardTypes;
}

double Config::getReactionForwardRate(unsigned i)
{
	return this->reactions[i]->forwardRate;
}

double Config::getReactionBackwardRate(unsigned i)
{
	return this->reactions[i]->backwardRate;
}