//
// Created by max on 8/8/16.
//

#include <gtest/gtest.h>
#include <random.hpp>
#include <folders/nussinov_folder.hpp>

#include <folders/brute_folder.hpp>

using namespace std;

TEST(NussinovFolder, BasicGGGGCCCC) {
	librnary::PrimeStructure struc = librnary::StringToPrimary("GGGGCCCC");
	librnary::NussinovFolder folder;
	int max = folder.Fold(struc.size(), [&](int i, int j) -> int {
		return librnary::ValidPair(struc[i], struc[j]) ? 1 : -1;
	});
	EXPECT_EQ(max, 4);
	auto match = folder.Traceback();
	EXPECT_EQ(librnary::MatchingToDotBracket(match), "(((())))");
}

TEST(NussinovFolder, BasicAUAGC) {
	librnary::PrimeStructure struc = librnary::StringToPrimary("AUAGC");
	librnary::NussinovFolder folder;
	int max = folder.Fold(struc.size(), [&](int i, int j) -> int {
		return librnary::ValidPair(struc[i], struc[j]) ? 1 : -1;
	});
	EXPECT_EQ(max, 2);
}

/*
 * A very simple scorer that returns the negation (for MFE) of the number of base pairs.
 */
struct SimpleScorer {
	librnary::PrimeStructure prim;
	int Score(const librnary::Surface &surf) const {
		assert(surf.IsExternalLoop() || librnary::ValidPair(prim[surf.PairI()], prim[surf.PairJ()]));
		int ans = 0;
		for (const auto &ss : surf.Children()) {
			ans += 1 + Score(ss);
		}
		return ans;
	}
	int ScoreExterior(const librnary::Surface &surf) const {
		return -Score(surf);
	}
	void SetRNA(const librnary::PrimeStructure &_rna) {
		prim = _rna;
	}
	int MaxMFE() const {
		return std::numeric_limits<int>::max() / 3;
	}
};

TEST(NussinovFolder, FuzzAgainstBrute) {
	SimpleScorer ss;
	librnary::StructureEnumerator se(0);
	librnary::BruteFolder bf(se);
	librnary::NussinovFolder folder;
	auto re = librnary::RandomEngineForTests();
	const int CASES = 1000, MAXN = 15;
	for (int tc = 0; tc < CASES; ++tc) {
		int N = (int) re() % MAXN;
		// Check same "MFE".
		auto rna = librnary::RandomPrimary(re, unsigned(N));
		ss.SetRNA(rna);
		auto bf_res = bf.FoldMFE(ss, rna);
		int folder_mfe = folder.Fold(rna.size(), [&](int i, int j) -> int {
			return librnary::ValidPair(rna[i], rna[j]) ? 1 : -1;
		}), bf_mfe = -get<0>(bf_res[0]);
		EXPECT_EQ(bf_mfe, folder_mfe);
		// Check that the traceback was in there.
		auto match = folder.Traceback();
		bool found = false;
		for (const auto &t : bf_res) {
			if (get<1>(t) == match) {
				found = true;
				break;
			}
		}
		EXPECT_TRUE(found);
	}
}