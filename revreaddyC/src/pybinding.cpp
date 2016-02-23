#include <string>
#include <vector>
#include <boost/python.hpp>
#include <boost/python/numeric.hpp>
#include <boost/python/extract.hpp>
#include "World.h"
#include "Config.h"
#include "Simulation.h"
#include "logging.h"


class ConfigWrap {
public:
	Config * config;
	ConfigWrap() { this->config = new Config();	}
	~ConfigWrap() { delete this->config; }
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

	void addParticle(boost::python::numeric::array initPos, unsigned particleTypeId) {
		std::vector<double> pos = {0., 0., 0.};
		try {
			pos[0] = boost::python::extract<double>(initPos[0]);
			pos[1] = boost::python::extract<double>(initPos[1]);
			pos[2] = boost::python::extract<double>(initPos[2]);
		} catch (...) {
			LOG_ERROR("Exception in accessing boost::python::numeric::array.")
			LOG_INFO("Particle is not added.")
			return;
		}
		world->addParticle(pos, particleTypeId);
	}

	boost::python::numeric::array getPosition(unsigned long index) {
		std::vector<double> pos = world->getPosition(index);
		boost::python::numeric::array arr = boost::python::numeric::array(
			boost::python::make_tuple(pos[0], pos[1], pos[2])
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
};

using namespace boost::python;

BOOST_PYTHON_MODULE(simpy) {
	numeric::array::set_module_and_type("numpy", "ndarray");
	class_<WorldWrap>("World", init<>())
		.def("addParticle", &WorldWrap::addParticle)
		.def("getPosition", &WorldWrap::getPosition);
	class_<ConfigWrap>("Config", init<>());
	class_<SimulationWrap>("Simulation", init<WorldWrap*, ConfigWrap*, std::string>());
};