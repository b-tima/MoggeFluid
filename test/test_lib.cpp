#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

extern "C" {
#include "example_lib.h"
}

TEST_GROUP(example_group){void setup(){} void teardown(){}};

TEST(example_group, example_test) {
	CHECK_EQUAL(3, example_testable_function(2));
}
