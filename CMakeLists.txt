cmake_minimum_required(VERSION 2.8)
project(revreaddy)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

find_package(PythonInterp)
find_package(PythonLibs REQUIRED)
find_package(Boost COMPONENTS python)

include_directories(${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS} include /usr/include/python2.7)
link_libraries(${Boost_LIBRARIES} ${PYTHON_LIBRARIES} gsl gslcblas m hdf5 hdf5_hl)
#include_directories(include)
#link_libraries(gsl gslcblas m hdf5 hdf5_hl)
# these are the C++ standalone sources
set(SOURCES
    src/Acceptance.cpp
    src/BinaryFile.cpp
    src/Config.cpp
    src/Conversion.cpp
    src/DoubleWellZ.cpp
    src/Energy.cpp
    src/Fusion3.cpp
    src/Geometry.cpp
    src/LennardJones.cpp
    src/MeanSquaredDisplacement.cpp
    src/Observable.cpp
    src/Particle.cpp
    src/ParticleInteraction.cpp
    src/ParticleNumbers.cpp
    src/ParticleType.cpp
    src/ProbabilityDensity.cpp
    src/RadialDistribution.cpp
    src/Random.cpp
    src/Reaction.cpp
    src/ReactionEvent.cpp
    src/Simulation.cpp
    src/SoftRepulsion.cpp
    src/Trajectory.cpp
    src/UnimolecularCandidate.cpp
    src/Wall.cpp
    src/World.cpp
)

python_add_module(simpy ${SOURCES} src/pybinding.cpp)
add_library(simcpp SHARED ${SOURCES})
add_library(simcppstandalone STATIC ${SOURCES})