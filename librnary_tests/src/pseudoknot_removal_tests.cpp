//
// Created by max on 8/9/16.
//

#include <gtest/gtest.h>
#include <secondary_structure.hpp>
#include <pseudoknot_removal.hpp>

using namespace std;

TEST(PseudoRemoval, SimplePseudo1) {
	librnary::Matching match = {2, 3, 0, 1};
	match = librnary::RemovePseudoknotsMaximizePairs(match);
	EXPECT_EQ(librnary::MatchingToBondPairs(match).size(), 1);
}

TEST(PseudoRemoval, SimplePseudo2) {
	librnary::Matching match = {0, 4, 5, 7, 1, 2, 6, 3};
	match = librnary::RemovePseudoknotsMaximizePairs(match);
	EXPECT_EQ(librnary::MatchingToBondPairs(match).size(), 1);
}

TEST(PseudoRemoval, SimplePseudo3) {
	librnary::Matching match = {0, 5, 4, 7, 1, 2, 6, 3};
	match = librnary::RemovePseudoknotsMaximizePairs(match);
	EXPECT_EQ(librnary::MatchingToBondPairs(match).size(), 2);
}

TEST(PseudoRemoval, CliqueOf3_1) {
	// .(([{))]}.
	librnary::Matching match = {0, 6, 5, 7, 8, 2, 1, 3, 4, 9};
	match = librnary::RemovePseudoknotsMaximizePairs(match);
	EXPECT_EQ(librnary::MatchingToBondPairs(match).size(), 2);
}

TEST(PseudoRemoval, NoPseudo) {
	librnary::Matching match = {1, 0};
	match = librnary::RemovePseudoknotsMaximizePairs(match);
	EXPECT_EQ(librnary::MatchingToBondPairs(match).size(), 1);
}
