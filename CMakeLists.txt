cmake_minimum_required(VERSION 2.8)
project(revreaddy)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

find_package(PythonInterp)
find_package(PythonLibs REQUIRED)
find_package(Boost COMPONENTS python)

include_directories(${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS} include /usr/include/python2.7)
link_libraries(${Boost_LIBRARIES} ${PYTHON_LIBRARIES} gsl gslcblas m hdf5 hdf5_hl)
file(GLOB SOURCES "src/*.cpp")

python_add_module(simPy ${SOURCES} sim.cpp)
#add_library(simC SHARED ${SOURCES})