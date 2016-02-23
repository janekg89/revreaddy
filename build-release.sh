mkdir -p  build/Release
cd build/Release
cmake ../.. \
	-G "Unix Makefiles" \
	-DCMAKE_BUILD_TYPE=Release \
	-DREVREADDY_LOGGING_SEVERITY:STRING=__ERROR__ \
	-DREVREADDY_RUN_TESTS:BOOL=OFF
make -j
cd ../..