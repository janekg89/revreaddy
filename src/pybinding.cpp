#include <string>
#include "World.h"
#include "Config.h"
#include "Simulation.h"

class ConfigWrap {
public:
	Config * config;
	ConfigWrap() { this->config = new Config();	}
	~ConfigWrap() { delete this->config; }
};

class WorldWrap {
public:
	World * world;
	WorldWrap() { this->world = new World(); }
	~WorldWrap() { delete this->world; }
};

class SimulationWrap {
public:
	Simulation * simulation;
	SimulationWrap(WorldWrap * worldWrap, ConfigWrap * configWrap, std::string whichImpl) {
		this->simulation = new Simulation(worldWrap->world, configWrap->config, whichImpl);
	}
	~SimulationWrap() { delete simulation; }
};

#include <boost/python.hpp>
using namespace boost::python;

BOOST_PYTHON_MODULE(simpy) {
	class_<WorldWrap>("World", init<>());
	class_<ConfigWrap>("Config", init<>());
	class_<SimulationWrap>("Simulation", init<WorldWrap*, ConfigWrap*, std::string>());
};