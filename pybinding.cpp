#include <string>
#include <vector>
#include <boost/python.hpp>
#include <boost/python/numeric.hpp>
#include <boost/python/extract.hpp>
#include "World.h"
#include "Config.h"
#include "Simulation.h"
#include "logging.h"

namespace bp = boost::python;

class ConfigWrap {
public:
	Config * config;
	ConfigWrap() { this->config = new Config();	}
	~ConfigWrap() { delete this->config; }
	double getTimestep() { return this->config->timestep; }
	void setTimestep(double timestep) { this->config->timestep = timestep; }
	double getKT() { return this->config->kT; }
	void setKT(double kT) {this->config->kT = kT; }
	bool getIsPeriodic() { return this->config->isPeriodic; }
	void setIsPeriodic(bool isPeriodic) { this->config->isPeriodic = isPeriodic; }
	double getBoxsize() {return this->config->boxsize; }
	void setBoxsize(double boxsize) { this->config->boxsize = boxsize; }

	void deleteAllParticleTypes() {
		this->config->deleteAllParticleTypes();
	}
	void new_Type(const std::string name, const double radius, const double diffusionConstant) {
		this->config->new_Type(name, radius, diffusionConstant);
	}
	unsigned getNumberOfParticleTypes() {
		return this->config->getNumberOfParticleTypes();
	}
	std::string getParticleTypeName(unsigned i) {
		return this->config->getParticleTypeName(i);
	}
	double getParticleTypeRadius(unsigned i) {
		return this->config->getParticleTypeRadius(i);
	}
	double getParticleTypeDiffusionConstant(unsigned i) {
		return this->config->getParticleTypeDiffusionConstant(i);
	}

	void deleteAllGeometries() {
		this->config->deleteAllGeometries();
	}
	void new_Wall(bp::numeric::array normal, bp::numeric::array point, double strength, bp::numeric::array particleTypeIds) {
		std::vector<double> norm = {0.,0.,0.};
		std::vector<double> pnt = {0.,0.,0.};
		std::vector<unsigned> pTypeIds;
		try {
			norm[0] = bp::extract<double>(normal[0]);
			norm[1] = bp::extract<double>(normal[1]);
			norm[2] = bp::extract<double>(normal[2]);
			pnt[0] = bp::extract<double>(point[0]);
			pnt[1] = bp::extract<double>(point[1]);
			pnt[2] = bp::extract<double>(point[2]);
			for (unsigned i=0; i<particleTypeIds.nelements(); ++i) {
				pTypeIds.push_back(bp::extract<unsigned>(particleTypeIds[i]));
			}
		} catch (...) {
			LOG_ERROR("Exception in accessing bp::numeric::array in new_Wall.")
			LOG_INFO("Wall geometry is not added.")
			return;
		}
		this->config->new_Wall(norm, pnt, strength, pTypeIds);
	}
	void new_DoubleWellZ(double distanceMinima,	double strength, bp::numeric::array particleTypeIds) {
		std::vector<unsigned> pTypeIds;
		try {
			for (unsigned i=0; i<particleTypeIds.nelements(); ++i) {
				pTypeIds.push_back(bp::extract<unsigned>(particleTypeIds[i]));
			}		
		} catch (...) {
			LOG_ERROR("Exception in accessing bp::numeric::array in new_DoubleWellZ.")
			LOG_INFO("DoubleWellZ geometry is not added.")
			return;
		}
		this->config->new_DoubleWellZ(distanceMinima, strength, pTypeIds);
	}

	void deleteAllInteractions() {
		this->config->deleteAllInteractions();
	}
	void new_SoftRepulsion(std::string name, bp::numeric::array affectedTuple, double repulsionStrength) {
		std::vector<unsigned> affTuple = {0,0};
		try {
			affTuple[0] = bp::extract<unsigned>(affectedTuple[0]);
			affTuple[1] = bp::extract<unsigned>(affectedTuple[1]);
		} catch (...) {
			LOG_ERROR("Exception in accessing bp::numeric::array in new_SoftRepulsion.")
			LOG_INFO("SoftRepulsion interaction is not added.")
			return;
		}
		this->config->new_SoftRepulsion(name, affTuple, repulsionStrength);
	}
	void new_LennardJones(std::string name,	bp::numeric::array affectedTuple, double epsilon) {
		std::vector<unsigned> affTuple = {0,0};
		try {
			affTuple[0] = bp::extract<unsigned>(affectedTuple[0]);
			affTuple[1] = bp::extract<unsigned>(affectedTuple[1]);
		} catch (...) {
			LOG_ERROR("Exception in accessing bp::numeric::array in new_LennardJones.")
			LOG_INFO("LennardJones interaction is not added.")
			return;
		}
		this->config->new_LennardJones(name, affTuple, epsilon);
	}
	unsigned getNumberInteractions() {
		return this->config->getNumberInteractions();
	}
	std::string getInteractionName(unsigned i) {
		return this->config->getInteractionName(i);
	}
	std::string getInteractionType(unsigned i) {
		return this->config->getInteractionType(i);
	}
	bp::numeric::array getInteractionAffectedTuple(unsigned i) {
		std::vector<unsigned> affTuple = this->config->getInteractionAffectedTuple(i);
		bp::numeric::array affectedTuple = bp::numeric::array(
			bp::make_tuple(affTuple[0], affTuple[1])
		);
		return affectedTuple;
	}
	bp::numeric::array getInteractionParameters(unsigned i) {
		std::vector<double> interactionParameters = this->config->getInteractionParameters(i);
		// create list and append the doubles to it, then convert the list to numeric array and return it
		bp::list tempList = bp::list();
		for (unsigned j=0; j<interactionParameters.size(); ++j) {
			tempList.append<double>(interactionParameters[j]);
		}
		bp::numeric::array interactionParams = bp::numeric::array(tempList); 
		return interactionParams;
	}
	double getInteractionCutoff(unsigned i) {
		return this->config->getInteractionCutoff(i);
	}

	void deleteAllReactions() {
		this->config->deleteAllReactions();
	}
	void new_Conversion(std::string name, unsigned forwardType, unsigned backwardType, double forwardRate, double backwardRate) {
		this->config->new_Conversion(name, forwardType, backwardType, forwardRate, backwardRate);
	}
	void new_Fusion(std::string name, unsigned forwardTypeA, unsigned forwardTypeB,	unsigned backwardTypeC,	double forwardRate,	double backwardRate, double reactionDistance) {
		this->config->new_Fusion(name, forwardTypeA, forwardTypeB, backwardTypeC, forwardRate, backwardRate, reactionDistance);
	}
	// TODO this is WIP as long as Fusion is configured manually
	void configureFusion(unsigned reactionIndex, bp::numeric::array interactionsIndices,	double inversePartition, double maxDistr, double radiiSum, double reactionRadiiSum,	double meanDistr, double inverseTemperature, double radiusA, double radiusB) {
		std::vector<unsigned> interactionsIndicesConverted;
		try {
			for (unsigned i=0; i<interactionsIndices.nelements(); ++i) {
				interactionsIndicesConverted.push_back(	bp::extract<unsigned>(interactionsIndices[i]) );
			}
		} catch (...) {
			LOG_ERROR("Exception in accessing bp::numeric::array in configureFusion.")
			LOG_INFO("Fusion is not configured.")
			return;			
		}
		this->config->configureFusion(reactionIndex, interactionsIndicesConverted, inversePartition, maxDistr, radiiSum, reactionRadiiSum, meanDistr, inverseTemperature, radiusA, radiusB);
	}
	unsigned getNumberReactions() {
		return this->config->getNumberReactions();
	}
	std::string getReactionName(unsigned i) {
		return this->config->getReactionName(i);
	}
	std::string getReactionType(unsigned i) {
		return this->config->getReactionType(i);
	}
	bp::numeric::array getReactionForwardTypes(unsigned i) {
		std::vector<unsigned> forwardTypes = this->config->getReactionForwardTypes(i);
		bp::list tempList = bp::list();
		for (unsigned j=0; j<forwardTypes.size(); ++j) {
			tempList.append<unsigned>(forwardTypes[j]);
		}
		bp::numeric::array forwTypes = bp::numeric::array(tempList);
		return forwTypes;
	}
	bp::numeric::array getReactionBackwardTypes(unsigned i) {
		std::vector<unsigned> backwardTypes = this->config->getReactionBackwardTypes(i);
		bp::list tempList = bp::list();
		for (unsigned j=0; j<backwardTypes.size(); ++j) {
			tempList.append<unsigned>(backwardTypes[j]);
		}
		bp::numeric::array backwTypes = bp::numeric::array(tempList);
		return backwTypes;
	}
	double getReactionForwardRate(unsigned i) {
		return this->config->getReactionForwardRate(i);
	}
	double getReactionBackwardRate(unsigned i) {
		return this->config->getReactionBackwardRate(i);
	}
};

class WorldWrap {
public:
	World * world;
	WorldWrap() {
		LOG_TRACE("Enter WorldWrap Constructor.")
		this->world = new World();
		LOG_TRACE("Leave WorldWrap Constructor.")
	}

	~WorldWrap() {
		LOG_TRACE("Enter WorldWrap Destructor.")
		delete this->world;
		LOG_TRACE("Leave WorldWrap Destructor.")
	}

	void addParticle(bp::numeric::array initPos, unsigned particleTypeId) {
		std::vector<double> pos = {0., 0., 0.};
		try {
			pos[0] = bp::extract<double>(initPos[0]);
			pos[1] = bp::extract<double>(initPos[1]);
			pos[2] = bp::extract<double>(initPos[2]);
		} catch (...) {
			LOG_ERROR("Exception in accessing bp::numeric::array in addParticle.")
			LOG_INFO("Particle is not added.")
			return;
		}
		world->addParticle(pos, particleTypeId);
	}

	bp::numeric::array getPosition(unsigned long index) {
		std::vector<double> pos = world->getPosition(index);
		bp::numeric::array arr = bp::numeric::array(
			bp::make_tuple(pos[0], pos[1], pos[2])
		);
		return arr;
	}
};

class SimulationWrap {
public:
	Simulation * simulation;
	SimulationWrap(WorldWrap * worldWrap, ConfigWrap * configWrap, std::string whichImpl) {
		LOG_TRACE("Enter SimulationWrap Constructor.")
		this->simulation = new Simulation(worldWrap->world, configWrap->config, whichImpl);
		LOG_TRACE("Leave SimulationWrap Constructor.")
	}
	~SimulationWrap() {
		LOG_TRACE("Enter SimulationWrap Destructor.")
		delete simulation;
		LOG_TRACE("Leave SimulationWrap Destructor.")
	}
	void run(const unsigned long maxTime) {
		this->simulation->run(maxTime);
	}
	bool getUseNeighborlist() {
		return this->simulation->getUseNeighborlist();
	}
	void setUseNeighborlist(const bool inUseNeighborlist) {
		this->simulation->setUseNeighborlist(inUseNeighborlist);
	}

	void writeAllObservablesToFile() {
		this->simulation->writeAllObservablesToFile();
	}
	std::string showObservables() {
		return this->simulation->showObservables();
	}
	void deleteAllObservables() {
		this->simulation->deleteAllObservables();
	}
	void new_Trajectory(unsigned long recPeriod, std::string filename) {
		this->simulation->new_Trajectory(recPeriod, filename);
	}
	void new_RadialDistribution(unsigned long recPeriod, std::string filename, bp::numeric::array ranges, bp::numeric::array considered) {
		std::vector<double> rangesConverted;
		std::vector< std::vector<unsigned> > consideredConverted;
		try {
			for (unsigned i=0; i<ranges.nelements(); ++i) {
				rangesConverted.push_back( bp::extract<double>(ranges[i]) );
			}			
			// we expect considered to be 2-dimensional with x elements in 
			// first dimension and two elements in second dimension
			for (unsigned i=0; i<considered.nelements(); ++i) {
				std::vector<unsigned> consideredTuple;
				consideredTuple.push_back( bp::extract<unsigned>(considered[i][0]) );
				consideredTuple.push_back( bp::extract<unsigned>(considered[i][1]) );
				consideredConverted.push_back(consideredTuple);
			}
		} catch (...) {
			LOG_ERROR("Exception in accessing bp::numeric::array new_RadialDistribution.")
			LOG_INFO("RadialDistribution geometry is not added.")
		}
		this->simulation->new_RadialDistribution(recPeriod, filename, rangesConverted, consideredConverted);
	}
	void new_MeanSquaredDisplacement(unsigned long recPeriod, std::string filename,	unsigned particleTypeId) {
		this->simulation->new_MeanSquaredDisplacement(recPeriod, filename, particleTypeId);
	}
	void new_ProbabilityDensity(unsigned long recPeriod, std::string filename, unsigned pTypeId, bp::numeric::array range, unsigned int coord) {
		std::vector<double> rangeConverted;
		try {
			for (unsigned i=0; i<range.nelements(); ++i) {
				rangeConverted.push_back( bp::extract<double>(range[i]) );
			}
		} catch (...) {
			LOG_ERROR("Exception in accessing bp::numeric::array new_ProbabilityDensity.")
			LOG_INFO("RadialDistribution geometry is not added.")
		}
		this->simulation->new_ProbabilityDensity(recPeriod, filename, pTypeId, rangeConverted, coord);
	}
	void new_Energy(unsigned long recPeriod, std::string filename) {
		this->simulation->new_Energy(recPeriod, filename);
	}
	void new_Acceptance(unsigned long recPeriod, std::string filename, bool reactionsOrDynamics) {
		this->simulation->new_Acceptance(recPeriod, filename, reactionsOrDynamics);
	}
	void new_ParticleNumbers(unsigned long recPeriod, std::string filename,	unsigned particleTypeId) {
		this->simulation->new_ParticleNumbers(recPeriod, filename, particleTypeId);
	}
};

using namespace boost::python;

BOOST_PYTHON_MODULE(revreaddyPy) {
	numeric::array::set_module_and_type("numpy", "ndarray");
	class_<WorldWrap>("World", init<>())
		.def("addParticle", &WorldWrap::addParticle)
		.def("getPosition", &WorldWrap::getPosition);
	class_<ConfigWrap>("Config", init<>())
		.def("getTimestep", &ConfigWrap::getTimestep)
		.def("setTimestep", &ConfigWrap::setTimestep)
		.def("getKT", &ConfigWrap::getKT)
		.def("setKT", &ConfigWrap::setKT)
		.def("getIsPeriodic", &ConfigWrap::getIsPeriodic)
		.def("setIsPeriodic", &ConfigWrap::setIsPeriodic)
		.def("getBoxsize", &ConfigWrap::getBoxsize)
		.def("setBoxsize", &ConfigWrap::setBoxsize)
		.def("deleteAllParticleTypes", &ConfigWrap::deleteAllParticleTypes)
		.def("new_Type", &ConfigWrap::new_Type)
		.def("getNumberOfParticleTypes", &ConfigWrap::getNumberOfParticleTypes)
		.def("getParticleTypeName", &ConfigWrap::getParticleTypeName)
		.def("getParticleTypeRadius", &ConfigWrap::getParticleTypeRadius)
		.def("getParticleTypeDiffusionConstant", &ConfigWrap::getParticleTypeDiffusionConstant)
		.def("deleteAllGeometries", &ConfigWrap::deleteAllGeometries)
		.def("new_Wall", &ConfigWrap::new_Wall)
		.def("new_DoubleWellZ", &ConfigWrap::new_DoubleWellZ)
		.def("new_SoftRepulsion", &ConfigWrap::new_SoftRepulsion)
		.def("new_LennardJones", &ConfigWrap::new_LennardJones)
		.def("getNumberInteractions", &ConfigWrap::getNumberInteractions)
		.def("getInteractionName", &ConfigWrap::getInteractionName)
		.def("getInteractionType", &ConfigWrap::getInteractionType)
		.def("getInteractionAffectedTuple",&ConfigWrap::getInteractionAffectedTuple)
		.def("getInteractionParameters", &ConfigWrap::getInteractionParameters)
		.def("getInteractionCutoff", &ConfigWrap::getInteractionCutoff)
		.def("deleteAllReactions", &ConfigWrap::deleteAllReactions)
		.def("new_Conversion", &ConfigWrap::new_Conversion)
		.def("new_Fusion", &ConfigWrap::new_Fusion)
		.def("configureFusion", &ConfigWrap::configureFusion)
		.def("getNumberReactions", &ConfigWrap::getNumberReactions)
		.def("getReactionName", &ConfigWrap::getReactionName)
		.def("getReactionType", &ConfigWrap::getReactionType)
		.def("getReactionForwardTypes", &ConfigWrap::getReactionForwardTypes)
		.def("getReactionBackwardTypes", &ConfigWrap::getReactionBackwardTypes)
		.def("getReactionForwardRate", &ConfigWrap::getReactionForwardRate)
		.def("getReactionBackwardRate", &ConfigWrap::getReactionBackwardRate);
	class_<SimulationWrap>("Simulation", init<WorldWrap*, ConfigWrap*, std::string>())
		.def("run", &SimulationWrap::run)
		.def("getUseNeighborlist", &SimulationWrap::getUseNeighborlist)
		.def("setUseNeighborlist", &SimulationWrap::setUseNeighborlist)
		.def("writeAllObservablesToFile", &SimulationWrap::writeAllObservablesToFile)
		.def("showObservables", &SimulationWrap::showObservables)
		.def("deleteAllObservables", &SimulationWrap::deleteAllObservables)
		.def("new_Trajectory", &SimulationWrap::new_Trajectory)
		.def("new_RadialDistribution", &SimulationWrap::new_RadialDistribution)
		.def("new_MeanSquaredDisplacement", &SimulationWrap::new_MeanSquaredDisplacement)
		.def("new_ProbabilityDensity", &SimulationWrap::new_ProbabilityDensity)
		.def("new_Energy", &SimulationWrap::new_Energy)
		.def("new_Acceptance", &SimulationWrap::new_Acceptance)
		.def("new_ParticleNumbers", &SimulationWrap::new_ParticleNumbers);
};