cmake_minimum_required(VERSION 2.8.12)
project(revreaddy)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

find_package(PythonInterp)
find_package(PythonLibs REQUIRED)
find_package(Boost COMPONENTS python log REQUIRED)

include_directories(
    ${Boost_INCLUDE_DIRS}
    ${PYTHON_INCLUDE_DIRS}
    /usr/include/python2.7
    include
    include/geometries
    include/interactions
    include/observables
    include/reactions
)
link_libraries(${Boost_LIBRARIES} ${PYTHON_LIBRARIES} gsl gslcblas m hdf5 hdf5_hl)
#include_directories(include)
#link_libraries(gsl gslcblas m hdf5 hdf5_hl)
# these are the C++ standalone sources
set(SOURCES
    src/geometries/Geometry.cpp
    src/geometries/DoubleWellZ.cpp
    src/geometries/Wall.cpp

    src/interactions/ParticleInteraction.cpp
    src/interactions/LennardJones.cpp
    src/interactions/SoftRepulsion.cpp

    src/observables/Observable.cpp
    src/observables/Acceptance.cpp
    src/observables/Energy.cpp
    src/observables/MeanSquaredDisplacement.cpp
    src/observables/ParticleNumbers.cpp
    src/observables/ProbabilityDensity.cpp
    src/observables/RadialDistribution.cpp
    src/observables/Trajectory.cpp

    src/reactions/Reaction.cpp
    src/reactions/Conversion.cpp
    src/reactions/Fusion.cpp

    src/BinaryFile.cpp
    src/Config.cpp
    src/Particle.cpp
    src/ParticleType.cpp
    src/Random.cpp
    src/ReactionEvent.cpp
    src/Simulation.cpp
    src/SimulationImpl.cpp
    src/World.cpp
)

if(CMAKE_BUILD_TYPE MATCHES Debug)
    add_definitions(-D__DEBUG__)
endif(CMAKE_BUILD_TYPE MATCHES Debug)

if(CMAKE_BUILD_TYPE MATCHES Release)
    add_definitions(-D__INFO__)
endif(CMAKE_BUILD_TYPE MATCHES Release)

add_definitions(-DBOOST_LOG_DYN_LINK)
python_add_module(simpy ${SOURCES} src/pybinding.cpp)
#add_library(simcpp SHARED ${SOURCES})
#add_library(simcppstandalone STATIC ${SOURCES})