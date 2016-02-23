mkdir -p build/Debug
cd build/Debug
cmake ../.. \
	-G "Unix Makefiles" \
	-DCMAKE_BUILD_TYPE=Debug \
	-DREVREADDY_LOGGING_SEVERITY:STRING=__DEBUG__ \
	-DREVREADDY_RUN_TESTS:BOOL=ON
make -j
cd ../..