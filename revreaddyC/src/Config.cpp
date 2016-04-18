/* Config.cpp */
#include "Config.h"

template<typename T, typename ...Args>
std::unique_ptr<T> make_unique( Args&& ...args ) {
	return std::unique_ptr<T>( new T( std::forward<Args>(args)... ) );
}

Config::Config() {
	this->timestep              = 0.001;
	this->kT                    = 1.;
	this->isPeriodic            = true;
	this->boxsize               = 10.;
}

Config::~Config() {
	this->deleteAllParticleTypes();
}

void Config::deleteAllParticleTypes() {
	// since interactions, reactions and geometries rely on certain particle types
	// it is safer to remove them here as well
	// otherwise we will have segfaults
	LOG_INFO("Deleting all particle types. This means deleting all geometries, reactions and interactions as well.")
	this->deleteAllGeometries();
	this->deleteAllReactions();
	this->deleteAllInteractions();
	this->particleTypes.clear();
}

void Config::new_Type(const std::string name, const double radius, const double diffusionConstant) {
	if ( radius < 0. ) { throw Exception("The particle radius must be non-negative."); }
	if ( diffusionConstant < 0. ) { throw Exception("The diffusionConstant must be non-negative"); }
	ParticleType pType(
		name,
		radius,
		diffusionConstant);
	this->particleTypes.push_back(pType);
}

unsigned Config::getNumberParticleTypes() { return this->particleTypes.size(); }

std::string Config::getParticleTypeName(unsigned i) {
	if (i >= this->particleTypes.size()) {
		throw Exception("Particle type does not exist.");
	}
	return this->particleTypes[i].name;
}

double Config::getParticleTypeRadius(unsigned i) {
	if (i >= this->particleTypes.size()) {
		throw Exception("Particle type does not exist.");
	}
	return this->particleTypes[i].radius;
}

double Config::getParticleTypeDiffusionConstant(unsigned i) {
	if (i >= this->particleTypes.size()) {
		throw Exception("Particle type does not exist.");
	}
	return this->particleTypes[i].diffusionConstant;
}

void Config::deleteAllGeometries() {
	/* Erase all the unique pointers, the geometries are thus deleted as well */
	LOG_INFO("Deleting all geometries.")
	this->geometries.clear();
}

void Config::new_Wall(std::string name, std::vector<double> normal, std::vector<double> point, double strength, std::vector<unsigned int> particleTypeIds) {
	// TODO test exceptions
	for (unsigned i=0; i<particleTypeIds.size(); ++i) {
		if (particleTypeIds[i] >= this->particleTypes.size()) {
			throw Exception("Particle type does not exist.");
		}
	}
	// sort the particle types so that binary search can be used during runtime
	std::sort(particleTypeIds.begin(), particleTypeIds.end());
	if (normal.size() != 3) {
		throw Exception("Dimension mismatch in normal vector.");
	}
	// maybe: abs(normal) - 1. < epsilon would be better
	if (1. != normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]) {
		throw Exception("Normal vector is not normalized.");
	}
	if (point.size() != 3) {
		throw Exception("Dimension mismatch in point vector.");
	}
	if (strength < 0.) {
		throw Exception("Repulsion strength must be non-negative.");
	}
	std::unique_ptr<Wall> wall = make_unique<Wall>(
		name,
		normal,
		point,
		strength,
		particleTypeIds);
	this->geometries.push_back( std::move(wall) );
	LOG_INFO("Wall geometry added to geometries.")
}

void Config::new_DoubleWellZ(std::string name, double distanceMinima, double strength, std::vector<unsigned int> particleTypeIds) {
	for (unsigned i=0; i<particleTypeIds.size(); ++i) {
		if (particleTypeIds[i] >= this->particleTypes.size()) {
			throw Exception("Particle type does not exist.");
		}
	}
	// sort the particle types so that binary search can be used during runtime
	std::sort(particleTypeIds.begin(), particleTypeIds.end());
	if (distanceMinima < 0.) {
		throw Exception("The distance of minima must be non-negative.");
	}
	if (strength < 0.) {
		throw Exception("The strength of the double well must be non-negative.");
	}
	std::unique_ptr<DoubleWellZ> well = make_unique<DoubleWellZ>(
		name,
		distanceMinima,
		strength,
		particleTypeIds);
	this->geometries.push_back( std::move(well) );
	LOG_INFO("DoubleWellZ geometry added to geometries.")
}

unsigned Config::getNumberGeometries() {
	return this->geometries.size();
}

std::string Config::getGeometryName(unsigned i) {
	if (i >= this->geometries.size()) {
		throw Exception("Geometry of given index does not exist.");
	}
	return this->geometries[i]->name;
}

std::string Config::getGeometryType(unsigned i) {
	if (i >= this->geometries.size()) {
		throw Exception("Geometry of given index does not exist.");
	}
	return this->geometries[i]->type;
}

std::vector<unsigned> Config::getGeometryAffected(unsigned i) {
	if (i >= this->geometries.size()) {
		throw Exception("Geometry of given index does not exist.");
	}
	return this->geometries[i]->particleTypeIds;
}

void Config::deleteAllInteractions() {
	/* Erase all the shared pointers, if some interactions are still 
	 * referenced by Reactions, the interaction is technically valid,
	 * but it will lead to unwanted behavior. Therefore delete all
	 * reactions here as well. */
	LOG_INFO("Deleting all interactions and reactions.")
	this->reactions.clear();
	this->interactions.clear();
}

void Config::new_SoftRepulsion(std::string name, std::array<unsigned,2> affectedTuple, double repulsionStrength) {
	LOG_TRACE("Enter Config new_SoftRepulsion")
	if ( ( affectedTuple[0] >= this->particleTypes.size() ) || (affectedTuple[1] >= this->particleTypes.size() ) ) {
		throw Exception("The given particle type(s) do not exist.");
	}
	if ( repulsionStrength <= 0. ) {
		throw Exception("The given repulsion strength is not positive.");
	}
	std::shared_ptr<SoftRepulsion> soft = std::make_shared<SoftRepulsion>(
		name,
		affectedTuple,
		repulsionStrength);
	// set cutoff
	soft->cutoff = this->particleTypes[affectedTuple[0]].radius + this->particleTypes[affectedTuple[1]].radius;
	this->interactions.push_back( std::move(soft) );
	LOG_INFO("SoftRepulsion interaction added to interactions")
}

void Config::new_LennardJones(std::string name, std::array<unsigned,2> affectedTuple, double epsilon) {
	if ( (affectedTuple[0] >= this->particleTypes.size() ) || (affectedTuple[1] >= this->particleTypes.size() ) ) {
		throw Exception("The given particle type(s) do not exist.");
	}
	if ( epsilon <= 0. ) { throw Exception("The given epsilon is not positive."); }
	std::shared_ptr<LennardJones> lj = std::make_shared<LennardJones>(
		name,
		affectedTuple,
		epsilon);
	// Set cutoff. This is hardcoded and is applied in implementation of LennardJones.
	lj->cutoff = 2.5 * ( this->particleTypes[affectedTuple[0]].radius + this->particleTypes[affectedTuple[1]].radius );
	this->interactions.push_back( std::move(lj) );
	LOG_INFO("LennardJones interaction added to interactions")
}

unsigned int Config::getNumberInteractions() {
	return this->interactions.size();
}

std::string Config::getInteractionName(unsigned i) {
	if (i >= this->interactions.size()) {
		throw Exception("Interaction does not exist.");
	}
	return this->interactions[i]->name;
}

std::string Config::getInteractionType(unsigned i) {
	if (i >= this->interactions.size()) {
		throw Exception("Interaction does not exist.");
	}
	return this->interactions[i]->type;
}

std::array<unsigned,2> Config::getInteractionAffectedTuple(unsigned i) {
	if (i >= this->interactions.size()) {
		throw Exception("Interaction does not exist.");
	}
	return this->interactions[i]->affectedTuple;
}

std::vector<double> Config::getInteractionParameters(unsigned i) {
	if (i >= this->interactions.size()) {
		throw Exception("Interaction does not exist.");
	}
	return this->interactions[i]->parameters;
}

double Config::getInteractionCutoff(unsigned i) {
	if (i >= this->interactions.size()) {
		throw Exception("Interaction does not exist.");
	}
	return this->interactions[i]->cutoff;
}

void Config::configureReactions() {
	for (unsigned i=0; i<this->reactions.size(); ++i) {
		if (this->reactions[i]->type == "Conversion") {
			/* Conversion needs no configuration */
		}
		else if (this->reactions[i]->type == "Fusion") {
			//this->configureFusion(i); // TODO
		}
	}
}

void Config::deleteAllReactions() {
	/* Erase all the unique pointers from reactions. Thus all
	 * reactions will be destroyed accordingly */
	LOG_INFO("Deleting all reactions.")
	this->reactions.clear();	
}

void Config::new_Conversion(
	std::string name,
	unsigned forwardType,
	unsigned backwardType,
	double forwardRate,
	double backwardRate)
{
	if (forwardType >= this->particleTypes.size()) {
		throw Exception("The forward particle type does not exist.");
	}
	if (backwardType >= this->particleTypes.size()) {
		throw Exception("The backward particle type does not exist.");
	}
	if (forwardRate < 0.) {
		throw Exception("The forwardRate must be non-negative.");
	}
	if (backwardRate < 0.) {
		throw Exception("The backwardRate must be non-negative.");
	}
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

void Config::new_Fusion(
	std::string name,
	unsigned forwardTypeA,
	unsigned forwardTypeB,
	unsigned backwardTypeC,
	double forwardRate,
	double backwardRate,
	double reactionDistance) 
{
	if (forwardTypeA >= this->particleTypes.size()) {
		throw Exception("The forward particle type A does not exist.");
	}
	if (forwardTypeB >= this->particleTypes.size()) {
		throw Exception("The forward particle type B does not exist.");
	}
	if (backwardTypeC >= this->particleTypes.size()) {
		throw Exception("The backward particle type C does not exist.");
	}
	if (forwardRate < 0.) {
		throw Exception("The forwardRate must be non-negative.");
	}
	if (backwardRate < 0.) {
		throw Exception("The backwardRate must be non-negative.");
	}
	if (reactionDistance < 0.) {
		throw Exception("The reaction distance must be non-negative.");
	}
	std::vector<unsigned> forwardTypes;
	forwardTypes.push_back(forwardTypeA);
	forwardTypes.push_back(forwardTypeB);
	std::vector<unsigned> backwardTypes;
	backwardTypes.push_back(backwardTypeC);
	std::unique_ptr<Fusion> fus = make_unique<Fusion>(
		name,
		forwardTypes,
		backwardTypes,
		forwardRate,
		backwardRate,
		reactionDistance);
	this->reactions.push_back( std::move(fus) );
}

/* That reactionIndex is of type Fusion has to be ensured by the caller */
void Config::configureFusion(
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
		 * forwarded to the Fusion. */
		consideredInteractions.push_back( this->interactions[ interactionsIndices[i] ]);
	}
	// TODO this still seems inappropriate
	Fusion * fus = dynamic_cast<Fusion*>( this->reactions[reactionIndex].get() );
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
	 * The reference was only borrowed to cast it to Fusion 
	 * and configure it properly */
}

unsigned Config::getNumberReactions() {
	return this->reactions.size();
}

std::string Config::getReactionName(unsigned i) {
	if (i >= this->reactions.size()) {
		throw Exception("Reaction does not exist.");
	}
	return this->reactions[i]->name;
}

std::string Config::getReactionType(unsigned i) {
	if (i >= this->reactions.size()) {
		throw Exception("Reaction does not exist.");
	}
	return this->reactions[i]->type;
}

std::vector<unsigned> Config::getReactionForwardTypes(unsigned i) {
	if (i >= this->reactions.size()) {
		throw Exception("Reaction does not exist.");
	}
	return this->reactions[i]->forwardTypes;
}

std::vector<unsigned> Config::getReactionBackwardTypes(unsigned i) {
	if (i >= this->reactions.size()) {
		throw Exception("Reaction does not exist.");
	}
	return this->reactions[i]->backwardTypes;
}

double Config::getReactionForwardRate(unsigned i) {
	if (i >= this->reactions.size()) {
		throw Exception("Reaction does not exist.");
	}
	return this->reactions[i]->forwardRate;
}

double Config::getReactionBackwardRate(unsigned i) {
	if (i >= this->reactions.size()) {
		throw Exception("Reaction does not exist.");
	}
	return this->reactions[i]->backwardRate;
}