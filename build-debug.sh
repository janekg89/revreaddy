mkdir bin/Debug
cd bin/Debug
cmake ../.. -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Debug
make -j
cd ../..