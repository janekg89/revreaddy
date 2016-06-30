#!/usr/bin/env bash
echo $PATH
mkdir -p build/Debug
cd build/Debug

CMAKE_FLAGS=""
CMAKE_FLAGS+=" -DCMAKE_BUILD_TYPE=Debug"
CMAKE_FLAGS+=" -DREVREADDY_LOGGING_SEVERITY:STRING=__DEBUG__"
CMAKE_FLAGS+=" -DREVREADDY_RUN_TESTS:BOOL=ON"
CMAKE_FLAGS+=" -DBOOST_ROOT=/home/mi/janekg89/miniconda2/envs/revreaddy"
CMAKE_FLAGS+=" -DPYTHON_LIBRARY=/home/mi/janekg89/miniconda2/envs/revreaddy/lib/libpython2.7.so"
CMAKE_FLAGS+=" -DPYTHON_INCLUDE_DIR=/home/mi/janekg89/miniconda2/envs/revreaddy/include/python2.7"
CMAKE_FLAGS+=" -DPYTHON_EXECUTABLE=/home/mi/janekg89/miniconda2/envs/revreaddy/bin/python2.7"
CMAKE_FLAGS+=" -DZLIB_ROOT=/home/mi/janekg89/miniconda2/envs/revreaddy"
CMAKE_FLAGS+=" -DCMAKE_RUNTIME_OUTPUT_DIRECTORY:STRING=."
export HDF5_ROOT="/home/mi/janekg89/miniconda2/envs/revreaddy"

cmake ../.. $CMAKE_FLAGS
make -j4 #VERBOSE=1
make test -j4 #VERBOSE=1
#ctest -VV
cd ../..

