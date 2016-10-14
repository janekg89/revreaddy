# find pthread lib
find_package(Threads REQUIRED)

# gather sources
file(GLOB TEST_SOURCES *.cpp)

link_directories(${GOOGLETEST_LIBRARY} ${GOOGLEMOCK_LIBRARY})

# message(FATAL_ERROR "BLALBLALBA ${HDF5_LIBRARIES}")

# set target and link needed libraries
add_executable(runTests ${TEST_SOURCES})
# add_dependencies(runTests ${REVREADDY_C_LIBRARY_NAME})
target_include_directories(runTests PUBLIC ${GOOGLETEST_INCLUDE} ${GOOGLEMOCK_INCLUDE} ${REVREADDY_C_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})
target_link_libraries(runTests gtest gtest_main gmock gmock_main ${REVREADDY_C_LIBRARY_NAME} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES}
        ${CMAKE_THREAD_LIBS_INIT} ${Boost_woPython_LIBRARIES} )

# add tests to the building process
add_test(NAME revreaddy-tests COMMAND runTests)