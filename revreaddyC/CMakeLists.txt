# these are the C++ standalone sources
file(GLOB REVREADDY_C_SOURCES
    src/*.cpp 
    src/geometries/*.cpp
    src/interactions/*.cpp
    src/observables/*.cpp
    src/reactions/*.cpp 
)

add_library(${REVREADDY_C_LIBRARY_NAME} SHARED ${REVREADDY_C_SOURCES})
target_include_directories(${REVREADDY_C_LIBRARY_NAME} PUBLIC
    ${Boost_INCLUDE_DIRS}
    ${PROJECT_SOURCE_DIR}/revreaddyC/include
    ${PROJECT_SOURCE_DIR}/revreaddyC/include/geometries
    ${PROJECT_SOURCE_DIR}/revreaddyC/include/interactions
    ${PROJECT_SOURCE_DIR}/revreaddyC/include/observables
    ${PROJECT_SOURCE_DIR}/revreaddyC/include/reactions
)
target_compile_definitions(${REVREADDY_C_LIBRARY_NAME} PUBLIC -DBOOST_LOG_DYN_LINK)
target_link_libraries(${REVREADDY_C_LIBRARY_NAME} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} gsl gslcblas m hdf5 hdf5_hl)
set_property(TARGET ${REVREADDY_C_LIBRARY_NAME} PROPERTY COMPILE_FLAGS -std=c++11)

# set logging defines
set(REVREADDY_LOGGING_SEVERITY __INFO__ CACHE STRING "The logging level")
if (REVREADDY_LOGGING_SEVERITY MATCHES __TRACE__)
    target_compile_definitions(${REVREADDY_C_LIBRARY_NAME} PUBLIC -D__TRACE__)
elseif (REVREADDY_LOGGING_SEVERITY MATCHES __DEBUG__)
    target_compile_definitions(${REVREADDY_C_LIBRARY_NAME} PUBLIC -D__DEBUG__)
elseif (REVREADDY_LOGGING_SEVERITY MATCHES __INFO__)
    target_compile_definitions(${REVREADDY_C_LIBRARY_NAME} PUBLIC -D__INFO__)
elseif (REVREADDY_LOGGING_SEVERITY MATCHES __WARNING__)
    target_compile_definitions(${REVREADDY_C_LIBRARY_NAME} PUBLIC -D__WARNING__)
endif (REVREADDY_LOGGING_SEVERITY MATCHES __TRACE__)