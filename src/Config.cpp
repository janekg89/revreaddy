/* Config.cpp */

#include "Config.h"
#define print(x) std::cout << x << std::endl;

template<typename T, typename ...Args>
std::unique_ptr<T> make_unique( Args&& ...args )
{
	return std::unique_ptr<T>( new T( std::forward<Args>(args)... ) );
}

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
	this->useNeighborList       = false;
	this->numberBoxes           = 1;
	this->reactionPropagation   = 0;
}

Config::~Config()
{
	this->deleteAllObservables();
	this->deleteAllGeometries();
	this->deleteAllReactions();
	this->deleteAllForces();
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
	for (unsigned i=0; i<this->observables.size(); ++i) {
		this->observables[i]->writeBufferToFile();
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
		for (unsigned i=0; i<this->observables.size(); ++i) {
			content += std::string (typeid(* this->observables[i]).name()) + " ";
		}
	}
	else {content += "empty";}
	content += "\n";
	return content;
}

void Config::deleteAllObservables()
{
	/* Erase all the unique pointers, the observables are thus deleted as well */
	this->observables.clear();
}

void Config::deleteLastObservable()
{
	/* delete the smart pointer, upon which the observable will be destroyed as well */
	if (this->observables.size() > 0) {
		this->observables.pop_back();
	}
}

void Config::new_Trajectory(unsigned long recPeriod, std::string filename)
{
	std::unique_ptr<Trajectory> obs = make_unique<Trajectory>(
		recPeriod,
		0,
		filename);
	this->observables.push_back( std::move(obs) );
}

void Config::new_RadialDistribution(
	unsigned long recPeriod,
	std::string filename,
	std::vector<double> ranges,
	std::vector< std::vector<unsigned> > considered)
{
	std::unique_ptr<RadialDistribution> rad = make_unique<RadialDistribution>(
		recPeriod,
		0,
		ranges,
		this->isPeriodic,
		this->boxsize,
		considered,
		filename);
	this->observables.push_back( std::move(rad) );
}

void Config::new_MeanSquaredDisplacement(
	unsigned long recPeriod,
	std::string filename,
	unsigned particleTypeId)
{
	std::unique_ptr<MeanSquaredDisplacement> msd = make_unique<MeanSquaredDisplacement>(
		recPeriod,
		0,
		world->activeParticles,
		particleTypeId,
		world->cumulativeRuntime,
		this->boxsize,
		filename);
	this->observables.push_back( std::move(msd) );
}

void Config::new_ProbabilityDensity(
	unsigned long recPeriod,
	std::string filename,
	unsigned particleTypeId,
	std::vector<double> range,
	unsigned coord)
{
	std::unique_ptr<ProbabilityDensity> prob = make_unique<ProbabilityDensity>(
		recPeriod,
		0,
		filename,
		particleTypeId,
		range,
		coord,
		this->world);
	this->observables.push_back( std::move(prob) );
}

void Config::new_Energy(unsigned long recPeriod, std::string filename)
{
	std::unique_ptr<Energy> ener = make_unique<Energy>(
		recPeriod,
		0,
		filename);
	this->observables.push_back( std::move(ener) );
}

void Config::new_Acceptance(
	unsigned long recPeriod,
	std::string filename,
	bool reactionsOrDynamics)
{
	std::unique_ptr<Acceptance> acc = make_unique<Acceptance>(
		recPeriod,
		0,
		filename,
		reactionsOrDynamics);
	this->observables.push_back( std::move(acc) );
}

void Config::new_ParticleNumbers(
	unsigned long recPeriod,
	std::string filename,
	unsigned particleTypeId)
{
	std::unique_ptr<ParticleNumbers> par = make_unique<ParticleNumbers>(
		recPeriod,
		0,
		filename,
		particleTypeId);
	this->observables.push_back( std::move(par) );
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

void Config::deleteAllForces()
{
	/* Erase all the shared pointers, if some interactions are still 
	 * referenced by Reactions, the interaction is technically valid,
	 * but it will lead to unwanted behavior. Therefore delete all
	 * reactions here as well. */
	this->possibleReactions.clear();
	this->possibleInteractions.clear();
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
	std::shared_ptr<SoftRepulsion> soft = std::make_shared<SoftRepulsion>(
		name,
		affectedTuple,
		repulsionStrength);
	// set cutoff correctly
	soft->cutoff = this->typeDict[affectedTuple[0]].radius 
	             + this->typeDict[affectedTuple[1]].radius;
	this->possibleInteractions.push_back( std::move(soft) );
	std::cout << "Info: SoftRepulsion interaction added to possibleInteractions"
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
	std::shared_ptr<LennardJones> lj = std::make_shared<LennardJones>(
		name,
		affectedTuple,
		epsilon);
	// set cutoff correctly
	lj->cutoff = 2.5 * ( this->typeDict[affectedTuple[0]].radius
	                   + this->typeDict[affectedTuple[1]].radius );
	this->possibleInteractions.push_back( std::move(lj) );
	std::cout << "Info: LennardJones interaction added to possibleInteractions"
	          << std::endl;
}

unsigned int Config::getNumberForces()
{
	return this->possibleInteractions.size();
}

std::string Config::getForceName(unsigned i)
{
	return this->possibleInteractions[i]->name;
}

std::string Config::getForceType(unsigned i)
{
	return this->possibleInteractions[i]->type;
}

std::vector<unsigned int> Config::getForceAffectedTuple(unsigned i)
{
	return this->possibleInteractions[i]->affectedTuple;
}

std::vector<double> Config::getForceParameters(unsigned i)
{
	return this->possibleInteractions[i]->parameters;
}

double Config::getForceCutoff(unsigned i)
{
	return this->possibleInteractions[i]->cutoff;
}

void Config::deleteAllReactions()
{
	/* Erase all the unique pointers from possibleReactions. Thus all
	 * reactions will be destroyed accordingly */
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
	std::unique_ptr<Conversion> conv = make_unique<Conversion>(
		name,
		forwardTypes,
		backwardTypes,
		forwardRate,
		backwardRate,
		this->random);
	this->possibleReactions.push_back( std::move(conv) );
}

void Config::new_Fusion(
	std::string name,
	unsigned forwardTypeA,
	unsigned forwardTypeB,
	unsigned backwardTypeC,
	double forwardRate,
	double backwardRate)
{
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
		this->random);
	this->possibleReactions.push_back( std::move(fus) );
}

/* That reactionIndex is of type Fusion has to be ensured by the caller */
void Config::configure_Fusion(
	unsigned reactionIndex,
	std::vector<unsigned> interactionsIndices,
	double inversePartition,
	double maxDistr,
	double radiiSum,
	double reactionRadiiSum,
	double meanDistr,
	double inverseTemperature)
{
	std::vector< std::shared_ptr<ParticleInteraction> > interactions;
	for (unsigned i=0; i<interactionsIndices.size(); i++) {
		/* push_back constructs a new shared ptr to be
		 * forwarded to the fusion. */
		interactions.push_back(
			this->possibleInteractions[ interactionsIndices[i] ]
		);
	}
	// TODO this still seems inappropriate
	Fusion * fus = dynamic_cast<Fusion*>( this->possibleReactions[reactionIndex].get() );
	fus->configure(
		interactions,
		inversePartition,
		maxDistr,
		radiiSum,
		reactionRadiiSum,
		meanDistr,
		inverseTemperature);	
	/* Don't delete fus, since it is uniquely owned by Config. 
	 * The reference was only borrowed to cast it to Fusion 
	 * and configure it properly */
}

void Config::new_Fusion2(
	std::string name,
	unsigned forwardTypeA,
	unsigned forwardTypeB,
	unsigned backwardTypeC,
	double forwardRate,
	double backwardRate)
{
	std::vector<unsigned> forwardTypes;
	forwardTypes.push_back(forwardTypeA);
	forwardTypes.push_back(forwardTypeB);
	std::vector<unsigned> backwardTypes;
	backwardTypes.push_back(backwardTypeC);
	std::unique_ptr<Fusion2> fus = make_unique<Fusion2>(
		name,
		forwardTypes,
		backwardTypes,
		forwardRate,
		backwardRate,
		this->random);
	this->possibleReactions.push_back( std::move(fus) );
}

/* That reactionIndex is of type Fusion2 has to be ensured by the caller */
void Config::configure_Fusion2(
	unsigned reactionIndex,
	std::vector<unsigned> interactionsIndices,
	double inversePartition,
	double maxDistr,
	double radiiSum,
	double reactionRadiiSum,
	double meanDistr,
	double inverseTemperature)
{
	std::vector< std::shared_ptr<ParticleInteraction> > interactions;
	for (unsigned i=0; i<interactionsIndices.size(); i++) {
		/* push_back constructs a new shared ptr to be
		 * forwarded to the Fusion2. */
		interactions.push_back(
			this->possibleInteractions[ interactionsIndices[i] ]
		);
	}
	// TODO this still seems inappropriate
	Fusion2 * fus = dynamic_cast<Fusion2*>( this->possibleReactions[reactionIndex].get() );
	fus->configure(
		interactions,
		inversePartition,
		maxDistr,
		radiiSum,
		reactionRadiiSum,
		meanDistr,
		inverseTemperature);	
	/* Don't delete fus, since it is uniquely owned by Config. 
	 * The reference was only borrowed to cast it to Fusion2 
	 * and configure it properly */
}

void Config::new_Fusion3(
	std::string name,
	unsigned forwardTypeA,
	unsigned forwardTypeB,
	unsigned backwardTypeC,
	double forwardRate,
	double backwardRate)
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
		this->random);
	this->possibleReactions.push_back( std::move(fus) );
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
	std::vector< std::shared_ptr<ParticleInteraction> > interactions;
	for (unsigned i=0; i<interactionsIndices.size(); i++) {
		/* push_back constructs a new shared ptr to be
		 * forwarded to the Fusion3. */
		interactions.push_back(
			this->possibleInteractions[ interactionsIndices[i] ]
		);
	}
	// TODO this still seems inappropriate
	Fusion3 * fus = dynamic_cast<Fusion3*>( this->possibleReactions[reactionIndex].get() );
	fus->configure(
		interactions,
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
	return this->possibleReactions.size();
}

std::string Config::getReactionName(unsigned i)
{
	return this->possibleReactions[i]->name;
}

std::string Config::getReactionType(unsigned i)
{
	return this->possibleReactions[i]->type;
}

std::vector<unsigned> Config::getReactionForwardTypes(unsigned i)
{
	return this->possibleReactions[i]->forwardTypes;
}

std::vector<unsigned> Config::getReactionBackwardTypes(unsigned i)
{
	return this->possibleReactions[i]->backwardTypes;
}

double Config::getReactionForwardRate(unsigned i)
{
	return this->possibleReactions[i]->forwardRate;
}

double Config::getReactionBackwardRate(unsigned i)
{
	return this->possibleReactions[i]->backwardRate;
}