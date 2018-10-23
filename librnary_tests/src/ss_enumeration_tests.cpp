//
// Created by max on 5/30/16.
//

#include <gtest/gtest.h>

#include "secondary_structure.hpp"
#include "ss_enumeration.hpp"
#include "random.hpp"

using namespace std;

set<librnary::Matching> ExtractMatchingSet(const vector<string> &strings) {
	set<librnary::Matching> ret;
	for (const auto &s : strings) {
		ret.insert(librnary::DotBracketToMatching(s));
	}
	return ret;
}

TEST(SSEnumeration, AAAUUU) {
	auto prim = librnary::StringToPrimary("AAAUUU");
	auto forward = librnary::ForwardBondTable(prim);
	vector<vector<int>> correct_forward = {{3, 4, 5}, {3, 4, 5}, {3, 4, 5}, {}, {}, {}};
	EXPECT_EQ(forward, correct_forward);
	auto backward = librnary::BackwardBondTable(prim);
	vector<vector<int>> correct_backward = {{}, {}, {}, {0, 1, 2}, {0, 1, 2}, {0, 1, 2}};
	EXPECT_EQ(backward, correct_backward);
	librnary::StructureEnumerator se(3);
	set<librnary::Matching> pairs_set;
	auto f = [&pairs_set](const librnary::Matching &pairs) {
		pairs_set.insert(pairs);
	};
	se.Enumerate<decltype(f)>(prim, f);
	EXPECT_EQ(ExtractMatchingSet({"(...).", "(....)", ".(...)", "......"}), pairs_set);
}

TEST(SSEnumeration, GGGGG) {
	auto prim = librnary::StringToPrimary("GGGGG");
	auto forward = librnary::ForwardBondTable(prim);
	vector<vector<int>> correct_forward = {{}, {}, {}, {}, {}};
	EXPECT_EQ(forward, correct_forward);
	auto backward = librnary::BackwardBondTable(prim);
	vector<vector<int>> correct_backward = {{}, {}, {}, {}, {}};
	EXPECT_EQ(backward, correct_backward);
	librnary::StructureEnumerator se(3);
	set<librnary::Matching> pairs_set;
	auto f = [&pairs_set](const librnary::Matching &pairs) {
		pairs_set.insert(pairs);
	};
	se.Enumerate<decltype(f)>(prim, f);
	EXPECT_EQ(ExtractMatchingSet({"....."}), pairs_set);
}

TEST(SSEnumeration, GCGCGC) {
	auto prim = librnary::StringToPrimary("GCGCGC");
	auto forward = librnary::ForwardBondTable(prim);
	vector<vector<int>> correct_forward = {{1, 3, 5}, {2, 4}, {3, 5}, {4}, {5}, {}};
	EXPECT_EQ(forward, correct_forward);
	auto backward = librnary::BackwardBondTable(prim);
	vector<vector<int>> correct_backward = {{}, {0}, {1}, {0, 2}, {1, 3}, {0, 2, 4}};
	EXPECT_EQ(backward, correct_backward);
	librnary::StructureEnumerator se(3);
	set<librnary::Matching> pairs_set;
	auto f = [&pairs_set](const librnary::Matching &pairs) {
		pairs_set.insert(pairs);
	};
	se.Enumerate<decltype(f)>(prim, f);
	EXPECT_EQ(ExtractMatchingSet({"......", "(....)"}), pairs_set);
}

TEST(SSEnumeration, FuzzContainsRandomValidSturctures) {
	const int CASES = 500, RNALEN = 10, TRIALS = 10, TESTS = 20;
	auto re = librnary::RandomEngineForTests();
	for (int tc = 0; tc < CASES; ++tc) {
		auto primary = librnary::RandomPrimary(re, RNALEN);
		unsigned trials = min((unsigned) librnary::BondPairs(primary).size(), (unsigned) TRIALS);
		set<librnary::Matching> ss_set;
		librnary::StructureEnumerator se(0);
		auto f = [&ss_set](const librnary::Matching &pairs) {
			ss_set.insert(pairs);
		};
		se.Enumerate<decltype(f)>(primary, f);
		for (int t = 0; t < TESTS; ++t) {
			auto ss = librnary::RandomMatching(primary, re, trials);
			EXPECT_TRUE(ss_set.count(ss));
		}
	}
}

TEST(SSEnumeration, FuzzLacksInvalidSturctures) {
	const int CASES = 500, RNALEN = 10, TRIALS = 10, TESTS = 20;
	auto re = librnary::RandomEngineForTests();
	for (int tc = 0; tc < CASES; ++tc) {
		auto primary = librnary::RandomPrimary(re, RNALEN);
		unsigned trials = min((unsigned) librnary::BondPairs(primary).size(), (unsigned) TRIALS);
		vector<librnary::BondPair> invalid;
		for (int i = 0; i < RNALEN; ++i) {
			for (int j = i + 1; j < RNALEN; ++j) {
				if (!librnary::ValidPair(primary[i], primary[j])) {
					invalid.push_back(librnary::BondPair(i, j));
				}
			}
		}
		set<librnary::Matching> ss_set;
		librnary::StructureEnumerator se(0);
		auto f = [&ss_set](const librnary::Matching &pairs) {
			ss_set.insert(pairs);
		};
		se.Enumerate<decltype(f)>(primary, f);
		for (int t = 0; t < TESTS; ++t) {
			auto ss = librnary::RandomMatching(primary, re, trials);
			// Add random invalid pair.
			auto bp = invalid[re() % invalid.size()];
			ss[bp.i] = bp.j;
			ss[bp.j] = bp.i;
			EXPECT_FALSE(ss_set.count(ss));
		}
	}
}