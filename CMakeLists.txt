cmake_minimum_required(VERSION 2.8.12)
project(revreaddy)
set(REVREADDY_C_LIBRARY_NAME "revreaddyC")
set(REVREADDY_PY_LIBRARY_NAME "revreaddyPy")

find_package(PythonInterp)
find_package(PythonLibs REQUIRED)
find_package(Boost COMPONENTS python log REQUIRED)



# target python extension
python_add_module(${REVREADDY_PY_LIBRARY_NAME} src/pybinding.cpp)
target_include_directories(${REVREADDY_PY_LIBRARY_NAME}
    ${Boost_INCLUDE_DIRS}
    ${PYTHON_INCLUDE_DIRS}
    revreaddy/include
    revreaddy/include/geometries
    revreaddy/include/interactions
    revreaddy/include/observables
    revreaddy/include/reactions
)
#target_compile_definitions(${REVREADDY_C_LIBRARY_NAME} -DBOOST_LOG_DYN_LINK)
#target_link_libraries(${REVREADDY_LIBRARY_C_NAME} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} gsl gslcblas m hdf5 hdf5_hl)
set_property (TARGET ${REVREADDY_PY_LIBRARY_NAME} PROPERTY COMPILE_FLAGS -std=c++11)

# this is where the c library is build 
add_subdirectory(revreaddy)

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
    add_dependencies(${REVREADDY_C_LIBRARY_NAME} googletest)
    enable_testing()
    include_directories(${GOOGLETEST_INCLUDE})
    add_subdirectory(tests)
else ()
    message(STATUS "The submodule googletest was not initialized or REVREADDY_RUN_TESTS was OFF, thus no testing.")
endif ()