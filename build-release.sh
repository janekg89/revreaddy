mkdir bin/Release
cd bin/Release
cmake ../.. -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release
make -j
cd ../..
