cmake_minimum_required(VERSION 2.8.12)
project(revreaddy)

# use c++11
#set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# build C library and Py wrapping independently
set(REVREADDY_C_LIBRARY_NAME "revreaddyC")
set(REVREADDY_PY_LIBRARY_NAME "revreaddyPy")
set(REVREADDY_SCRIPT_NAME "revreaddy")

# py extension is build here, so find required packages
find_package(PythonInterp REQUIRED)
find_package(PythonLibs REQUIRED)
find_package(Boost COMPONENTS python log filesystem system REQUIRED)
find_package(HDF5 COMPONENTS C CXX HL)

# python toplayer
add_custom_target(${REVREADDY_SCRIPT_NAME} SOURCES revreaddy.py DEPENDS ${REVREADDY_PY_LIBRARY_NAME})
# target python extension
python_add_module(${REVREADDY_PY_LIBRARY_NAME} pybinding.cpp)
# c library must be built before py library can be built
add_dependencies(${REVREADDY_PY_LIBRARY_NAME} ${REVREADDY_C_LIBRARY_NAME})
set(REVREADDY_C_INCLUDE_DIRS
    ${PROJECT_SOURCE_DIR}/revreaddyC/include
    ${PROJECT_SOURCE_DIR}/revreaddyC/include/geometries
    ${PROJECT_SOURCE_DIR}/revreaddyC/include/interactions
    ${PROJECT_SOURCE_DIR}/revreaddyC/include/observables
    ${PROJECT_SOURCE_DIR}/revreaddyC/include/reactions
    ${PROJECT_SOURCE_DIR}/revreaddyC/include/simulationImpls
)

target_include_directories(${REVREADDY_PY_LIBRARY_NAME} PUBLIC
    ${Boost_INCLUDE_DIRS}
    ${PYTHON_INCLUDE_DIRS}
    ${REVREADDY_C_INCLUDE_DIRS}
)
target_link_libraries(${REVREADDY_PY_LIBRARY_NAME} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} ${REVREADDY_C_LIBRARY_NAME})

# this is where the c library is build 
add_subdirectory(revreaddyC)

# tests for the c library
set(REVREADDY_RUN_TESTS ON CACHE BOOL "Build and run tests.")
if (EXISTS "${CMAKE_SOURCE_DIR}/libraries/googletest/CMakeLists.txt" AND REVREADDY_RUN_TESTS)
    message(STATUS "The submodule googletest was initialized, thus enable testing")

    include(ExternalProject)

    # get the headers and libs of googletest
    set(GOOGLETEST_ROOT ${CMAKE_SOURCE_DIR}/libraries/googletest)
    set(GOOGLETEST_INCLUDE ${GOOGLETEST_ROOT}/googletest/include)
    set(GOOGLETEST_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/libraries/googletest)
    # if gtest should be build without gmock, the libraries will be in libraries/googletest/googletest
    set(GOOGLETEST_LIBRARY ${GOOGLETEST_BINARY_DIR}/googlemock/gtest)

    # get the headers and libs of googlemock (which is built together with googletest)
    set(GOOGLEMOCK_ROOT ${CMAKE_SOURCE_DIR}/libraries/googletest)
    set(GOOGLEMOCK_INCLUDE ${GOOGLEMOCK_ROOT}/googlemock/include)
    set(GOOGLEMOCK_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/libraries/googletest)
    set(GOOGLEMOCK_LIBRARY ${GOOGLEMOCK_BINARY_DIR}/googlemock)

    ExternalProject_Add(
        googletest
        DOWNLOAD_COMMAND ""
        SOURCE_DIR ${GOOGLETEST_ROOT}
        BINARY_DIR ${GOOGLETEST_BINARY_DIR}
        CMAKE_CACHE_ARGS "-DBUILD_GTEST:bool=ON" "-DBUILD_GMOCK:bool=ON"
        "-Dgtest_force_shared_crt:bool=ON"
        INSTALL_COMMAND ""
    )
    add_dependencies(${REVREADDY_C_LIBRARY_NAME} googletest)

    enable_testing()
    include_directories(${GOOGLETEST_INCLUDE} ${GOOGLEMOCK_INCLUDE})
    add_subdirectory(revreaddyC/tests)
    add_subdirectory(revreaddyC/tests/unittests)
else ()
    message(STATUS "The submodule googletest was not initialized or REVREADDY_RUN_TESTS was OFF, thus no testing.")
endif ()

# tests for the python wrapping
if(REVREADDY_RUN_TESTS)
    add_test(NAME python-wrapping COMMAND ${PYTHON_EXECUTABLE} test_revreaddy.py)
endif(REVREADDY_RUN_TESTS)

# when everything is built move the convenience layer script to the build directory
file(COPY revreaddy.py DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(COPY libraries/print_table.py DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
# same for tests
file(COPY test_revreaddy.py DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
