mkdir -p build/Debug
cd build/Debug
#PYTHON_EXECUTABLE=/usr/bin/python2.7 \
#PYTHON_LIBRARY=/usr/lib/python2.7/config-x86_64-linux-gnu/libpython2.7.so \
#PYTHON_INCLUDE=/usr/include/python2.7 \
cmake ../.. \
	-G "Unix Makefiles" \
	-DCMAKE_BUILD_TYPE=Debug \
	-DREVREADDY_LOGGING_SEVERITY:STRING=__DEBUG__ \
	-DREVREADDY_RUN_TESTS:BOOL=ON
make -j
make test -j
cd ../..