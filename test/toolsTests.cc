/* Generated file, do not edit */

#ifndef CXXTEST_RUNNING
#define CXXTEST_RUNNING
#endif

#define _CXXTEST_HAVE_STD
// MakeDepend: cflags CXXTEST
#include <cxxtest/TestListener.h>
#include <cxxtest/TestTracker.h>
#include <cxxtest/TestRunner.h>
#include <cxxtest/RealDescriptions.h>
#include <cxxtest/ErrorPrinter.h>

int main() {
 return CxxTest::ErrorPrinter().run();
}
#include "test/toolsTests.h"

static toolsTests suite_toolsTests;

static CxxTest::List Tests_toolsTests = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_toolsTests( "test/toolsTests.h", 10, "toolsTests", suite_toolsTests, Tests_toolsTests );

static class TestDescription_toolsTests_test_QualityTranslation : public CxxTest::RealTestDescription {
public:
 TestDescription_toolsTests_test_QualityTranslation() : CxxTest::RealTestDescription( Tests_toolsTests, suiteDescription_toolsTests, 38, "test_QualityTranslation" ) {}
 void runTest() { suite_toolsTests.test_QualityTranslation(); }
} testDescription_toolsTests_test_QualityTranslation;

static class TestDescription_toolsTests_test_SumLog2 : public CxxTest::RealTestDescription {
public:
 TestDescription_toolsTests_test_SumLog2() : CxxTest::RealTestDescription( Tests_toolsTests, suiteDescription_toolsTests, 57, "test_SumLog2" ) {}
 void runTest() { suite_toolsTests.test_SumLog2(); }
} testDescription_toolsTests_test_SumLog2;

static class TestDescription_toolsTests_testBinarySearchFunctional : public CxxTest::RealTestDescription {
public:
 TestDescription_toolsTests_testBinarySearchFunctional() : CxxTest::RealTestDescription( Tests_toolsTests, suiteDescription_toolsTests, 118, "testBinarySearchFunctional" ) {}
 void runTest() { suite_toolsTests.testBinarySearchFunctional(); }
} testDescription_toolsTests_testBinarySearchFunctional;

static class TestDescription_toolsTests_testLowerBoundNormLog2 : public CxxTest::RealTestDescription {
public:
 TestDescription_toolsTests_testLowerBoundNormLog2() : CxxTest::RealTestDescription( Tests_toolsTests, suiteDescription_toolsTests, 139, "testLowerBoundNormLog2" ) {}
 void runTest() { suite_toolsTests.testLowerBoundNormLog2(); }
} testDescription_toolsTests_testLowerBoundNormLog2;

#include <cxxtest/Root.cpp>
