mkdir -p  build/Release
cd build/Release
cmake ../.. -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DREVREADDY_LOGGING_SEVERITY:STRING=__ERROR__
make -j
cd ../..