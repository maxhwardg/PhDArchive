//
// Created by max on 5/5/16.
//

/*
 * Contains tests for basic secondary structure tools.
 */

#include <gtest/gtest.h>

#include "secondary_structure.hpp"
#include "ss_enumeration.hpp"
#include "random.hpp"

using namespace std;

namespace librnary {

// Tests if converting from dot bracket to Matching to BondPairs and back doesn't change the input.
TEST(SecondaryStructure, DBToPairsToBondsToPairsToDB) {
	vector<string> dbs = {"...", "(((...)))", "(...)(((...)))",
						  "..((((..(((..(...).)))..(.)...))))", "........(((((((()))))))).."};
	for (const string &db : dbs) {
		auto pairs = DotBracketToMatching(db);
		auto bonds = MatchingToBondPairs(pairs);
		auto p2 = BondPairsToMatching(bonds, static_cast<unsigned>(pairs.size()));
		EXPECT_EQ(MatchingToDotBracket(p2), db);
	}
}

// Fuzz tests if converting from dot bracket to BondPairs or dot bracket and back doesn't change the input.
TEST(SecondaryStructure, FuzzyPairsToBondsToPairsToDBToPairs) {
	const int TEST_CASES = 100, RNA_SIZE = 100, TRIALS = 500;
	default_random_engine re = RandomEngineForTests();
	for (int tc = 0; tc < TEST_CASES; ++tc) {
		auto prim = RandomPrimary(re, RNA_SIZE);
		auto pairs = RandomMatching(prim, re, TRIALS);
		EXPECT_EQ(DotBracketToMatching(MatchingToDotBracket(pairs)), pairs);
		EXPECT_EQ(BondPairsToMatching(MatchingToBondPairs(pairs), RNA_SIZE), pairs);
	}
}

// Tests if ValidPair works correctly using a truth table.
TEST(SecondaryStructure, ValidPair) {
	vector<vector<bool>> truth_table =
		{   //          A, U, G, C
			/* A */    {0, 1, 0, 0},
			/* U */    {1, 0, 1, 0},
			/* G */    {0, 1, 0, 1},
			/* C */    {0, 0, 1, 0}
		};
	ASSERT_EQ(truth_table.size(), NUMBASES);
	for (int i = 0; i < NUMBASES; ++i) {
		for (int j = 0; j < NUMBASES; ++j) {
			EXPECT_EQ(truth_table[i][j],
					  ValidPair(static_cast<Base>(i), static_cast<Base>(j)));
		}
	}
}

// Tests if WatsonCrick works correctly using a truth table.
TEST(SecondaryStructure, WatsonCrick) {
	vector<vector<bool>> truth_table =
		{   //          A, U, G, C
			/* A */    {0, 1, 0, 0},
			/* U */    {1, 0, 0, 0},
			/* G */    {0, 0, 0, 1},
			/* C */    {0, 0, 1, 0}
		};
	ASSERT_EQ(truth_table.size(), NUMBASES);
	for (int i = 0; i < NUMBASES; ++i) {
		for (int j = 0; j < NUMBASES; ++j) {
			EXPECT_EQ(truth_table[i][j],
					  WatsonCrick(static_cast<Base>(i), static_cast<Base>(j)));
		}
	}
}

// Tests if pseudoknots are detected.
TEST(SecondaryStructure, DetectPseudoknot) {
	Matching pairs = {2, 3, 0, 1};
	EXPECT_TRUE(ContainsPseudoknot(pairs));
	pairs = {5, 2, 1, 6, -1, 0, 3};
	EXPECT_TRUE(ContainsPseudoknot(pairs));
}

// Detects if no pseudoknots are found when there are none.
TEST(SecondaryStructure, DetectPseudoknotNoFalsePositives) {
	Matching pairs = {3, 2, 1, 0};
	EXPECT_FALSE(ContainsPseudoknot(pairs));
	pairs = {0, 1, 2, 3, 4};
	EXPECT_FALSE(ContainsPseudoknot(pairs));
	pairs = DotBracketToMatching("..(((..)..)..((..(.)...)))..(((...)))");
	EXPECT_FALSE(ContainsPseudoknot(pairs));
}

TEST(SecondaryStructure, MustBeLonelyPair) {
	auto rna = StringToPrimary("AUCGAUGCGUAGGGGC");
	EXPECT_TRUE(MustBeLonelyPair(rna, 4, 5, 3));
	EXPECT_FALSE(MustBeLonelyPair(rna, 3, 7, 3));
	EXPECT_TRUE(MustBeLonelyPair(rna, 1, 10, 3));
	EXPECT_TRUE(MustBeLonelyPair(rna, 7, 11, 3));
	EXPECT_FALSE(MustBeLonelyPair(rna, 2, 8, 3));
}

TEST(SecondaryStructure, ExtractStems) {
	auto match = DotBracketToMatching("((..((...))..(....).))..((((.((...))))))");
	auto stems = ExtractStems(match);
	multiset<int> sizes = {2,2,1,4,2};
	for (const auto &s : stems) {
		EXPECT_GT(sizes.count(s.num_pairs), 0);
		sizes.erase(sizes.find(s.num_pairs));
	}
}

}