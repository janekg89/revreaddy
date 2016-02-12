#include "World.h"

#include <boost/python.hpp>
using namespace boost::python;

BOOST_PYTHON_MODULE(simPy) {
	class_<World>("World")
		.def(init<>())
		.def("addParticle", &World::addParticle)
	;
};