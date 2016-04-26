echo $PATH
mkdir -p build/Release
cd build/Release

CMAKE_FLAGS=""
CMAKE_FLAGS+=" -DCMAKE_BUILD_TYPE=Release"
CMAKE_FLAGS+=" -DREVREADDY_LOGGING_SEVERITY:STRING=__ERROR__"
CMAKE_FLAGS+=" -DREVREADDY_RUN_TESTS:BOOL=OFF"
CMAKE_FLAGS+=" -DBOOST_ROOT=/home/chris/miniconda2/envs/revreaddy"
CMAKE_FLAGS+=" -DPYTHON_LIBRARY=/home/chris/miniconda2/envs/revreaddy/lib/libpython2.7.so"
CMAKE_FLAGS+=" -DPYTHON_INCLUDE_DIR=/home/chris/miniconda2/envs/revreaddy/include/python2.7"
CMAKE_FLAGS+=" -DPYTHON_EXECUTABLE=/home/chris/miniconda2/envs/revreaddy/bin/python2.7"
CMAKE_FLAGS+=" -DZLIB_ROOT=/home/chris/miniconda2/envs/revreaddy"
CMAKE_FLAGS+=" -DCMAKE_RUNTIME_OUTPUT_DIRECTORY:STRING=."

cmake ../.. $CMAKE_FLAGS
make -j
cd ../..
