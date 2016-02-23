cmake_minimum_required(VERSION 2.8.12)
project(revreaddy)

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

# set logging defines
set(REVREADDY_LOGGING_SEVERITY __INFO__ CACHE STRING "The logging level")
if (REVREADDY_LOGGING_SEVERITY MATCHES __TRACE__)
    add_definitions(-D__TRACE__)
elseif (REVREADDY_LOGGING_SEVERITY MATCHES __DEBUG__)
    add_definitions(-D__DEBUG__)
elseif (REVREADDY_LOGGING_SEVERITY MATCHES __INFO__)
    add_definitions(-D__INFO__)
elseif (REVREADDY_LOGGING_SEVERITY MATCHES __WARNING__)
    add_definitions(-D__WARNING__)
endif (REVREADDY_LOGGING_SEVERITY MATCHES __TRACE__)

# target python extension
add_definitions(-DBOOST_LOG_DYN_LINK)
python_add_module(revreaddyPyEx ${SOURCES} src/pybinding.cpp)
target_link_libraries(revreaddyPyEx ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} gsl gslcblas m hdf5 hdf5_hl)
set_property (TARGET revreaddyPyEx PROPERTY COMPILE_FLAGS -std=c++11)

# tests
set(REVREADDY_RUN_TESTS ON CACHE BOOL "Build and run tests.")
if (EXISTS "${CMAKE_SOURCE_DIR}/libraries/googletest/CMakeLists.txt" AND REVREADDY_RUN_TESTS)
    message(STATUS "The submodule googletest was initialized, thus enable testing")

    include(ExternalProject)

    set(GOOGLETEST_ROOT ${CMAKE_SOURCE_DIR}/libraries/googletest)
    set(GOOGLETEST_INCLUDE ${GOOGLETEST_ROOT}/googletest/include)
    set(GOOGLETEST_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/libraries/googletest)
    set(GOOGLETEST_LIBRARY ${GOOGLETEST_BINARY_DIR}/googletest)

    ExternalProject_Add(
            googletest
            DOWNLOAD_COMMAND ""
            SOURCE_DIR ${GOOGLETEST_ROOT}
            BINARY_DIR ${GOOGLETEST_BINARY_DIR}
            CMAKE_CACHE_ARGS "-DBUILD_GTEST:bool=ON" "-DBUILD_GMOCK:bool=OFF"
            "-Dgtest_force_shared_crt:bool=ON"
            INSTALL_COMMAND ""
    )
    add_dependencies(${READDY_LIBRARY_NAME} googletest)
    enable_testing()
    include_directories(${GOOGLETEST_INCLUDE})
    add_subdirectory(tests)
else ()
    message(STATUS "The submodule googletest was not initialized or REVREADDY_RUN_TESTS was OFF, thus no testing.")
endif ()