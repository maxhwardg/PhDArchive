//
// Created by max on 5/5/16.
//

#include "gtest/gtest.h"

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	int returnValue;
	// Do setup here.
	returnValue = RUN_ALL_TESTS();
	// Do teardown here.
	return returnValue;
}