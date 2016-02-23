

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

add_library(${REVREADDY_LIBRARY_NAME} STATIC SOURCES)
target_include_directories(${REVREADDY_LIBRARY_NAME}
    ${Boost_INCLUDE_DIRS}
    ${PYTHON_INCLUDE_DIRS}
    /usr/include/python2.7
    include
    include/geometries
    include/interactions
    include/observables
    include/reactions
)
target_compile_definitions(${REVREADDY_LIBRARY_NAME} -DBOOST_LOG_DYN_LINK)
target_link_libraries(${REVREADDY_LIBRARY_NAME} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} gsl gslcblas m hdf5 hdf5_hl)
set_property (TARGET ${REVREADDY_LIBRARY_NAME} PROPERTY COMPILE_FLAGS -std=c++11)

# set logging defines
set(REVREADDY_LOGGING_SEVERITY __INFO__ CACHE STRING "The logging level")
if (REVREADDY_LOGGING_SEVERITY MATCHES __TRACE__)
    target_compile_definitions(${REVREADDY_LIBRARY_NAME} -D__TRACE__)
elseif (REVREADDY_LOGGING_SEVERITY MATCHES __DEBUG__)
    target_compile_definitions(${REVREADDY_LIBRARY_NAME} -D__DEBUG__)
elseif (REVREADDY_LOGGING_SEVERITY MATCHES __INFO__)
    target_compile_definitions(${REVREADDY_LIBRARY_NAME} -D__INFO__)
elseif (REVREADDY_LOGGING_SEVERITY MATCHES __WARNING__)
    target_compile_definitions(${REVREADDY_LIBRARY_NAME} -D__WARNING__)
endif (REVREADDY_LOGGING_SEVERITY MATCHES __TRACE__)