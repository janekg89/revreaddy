/* Generated file, do not edit */

#ifndef CXXTEST_RUNNING
#define CXXTEST_RUNNING
#endif

#define _CXXTEST_HAVE_STD
#include <cxxtest/TestListener.h>
#include <cxxtest/TestTracker.h>
#include <cxxtest/TestRunner.h>
#include <cxxtest/RealDescriptions.h>
#include <cxxtest/TestMain.h>
#include <cxxtest/ErrorPrinter.h>

int main( int argc, char *argv[] ) {
 int status;
    CxxTest::ErrorPrinter tmp;
    CxxTest::RealWorldDescription::_worldName = "cxxtest";
    status = CxxTest::Main< CxxTest::ErrorPrinter >( tmp, argc, argv );
    return status;
}
bool suite_ParticleTest_init = false;
#include "/home/chris/workspace/revreaddy/tests/unittests/ParticleTest.h"

static ParticleTest suite_ParticleTest;

static CxxTest::List Tests_ParticleTest = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_ParticleTest( "tests/unittests/ParticleTest.h", 10, "ParticleTest", suite_ParticleTest, Tests_ParticleTest );

static class TestDescription_suite_ParticleTest_testMove : public CxxTest::RealTestDescription {
public:
 TestDescription_suite_ParticleTest_testMove() : CxxTest::RealTestDescription( Tests_ParticleTest, suiteDescription_ParticleTest, 13, "testMove" ) {}
 void runTest() { suite_ParticleTest.testMove(); }
} testDescription_suite_ParticleTest_testMove;

#include <cxxtest/Root.cpp>
const char* CxxTest::RealWorldDescription::_worldName = "cxxtest";
