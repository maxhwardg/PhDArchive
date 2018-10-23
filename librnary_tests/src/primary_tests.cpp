//
// Created by max on 5/5/16.
//

/*
 * Contains tests for the primary structure components of librnary.
 */

#include <gtest/gtest.h>

#include "primary_structure.hpp"
#include "random.hpp"

using namespace std;

/*
 * Tests if an RNA string converted to PrimeStructure
 * then back to a string is equal to the original string.
 */
TEST(PrimaryStructure, StringToPrimaryToString) {
	string rna = "AUGCUAUUGCUAGUC";
	auto prim = librnary::StringToPrimary(rna);
	EXPECT_EQ(librnary::PrimaryToString(prim), rna);
}

/*
 * Fuzz tests if a PrimeStructure converted to a string
 * then back to a PrimeStructure is equal to the original PrimeStructure.
 */
TEST(PrimaryStructure, PrimaryToStringToPrimary) {
	const int CASES = 100, LENGTH = 100;
	default_random_engine re = librnary::RandomEngineForTests();
	for (int t = 0; t < CASES; ++t) {
		auto prim = librnary::RandomPrimary(re, LENGTH);
		string rna = librnary::PrimaryToString(prim);
		EXPECT_EQ(librnary::StringToPrimary(rna), prim);
	}
}