# find pthread lib
find_package(Threads REQUIRED)

# gather unittest sources
file(GLOB UNITTEST_SOURCES *.cpp)

link_directories(${GOOGLETEST_LIBRARY} ${GOOGLEMOCK_LIBRARY})

# set target and link libraries
add_executable(runUnitTests ${UNITTEST_SOURCES})
add_dependencies(runUnitTests ${REVREADDY_C_LIBRARY_NAME})
target_include_directories(runUnitTests PUBLIC ${GOOGLETEST_INCLUDE} ${GOOGLEMOCK_INCLUDE} ${REVREADDY_C_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
target_link_libraries(runUnitTests gtest gtest_main gmock gmock_main ${REVREADDY_C_LIBRARY_NAME}
        ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})

# add tests so that "make -test" runs them
add_test(NAME revreaddy-unittests COMMAND runUnitTests)