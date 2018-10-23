//
// Created by max on 8/9/17.
//

#include <gtest/gtest.h>

#include <string>

#include "multi_loop.hpp"

using namespace std;

TEST(MultiLoop, ExtractLoopRegions1) {
	string db = "..((..(...)..(((...(((...))).)))..(...)))..(...)";
	auto match = librnary::DotBracketToMatching(db);
	librnary::SSTree sst(match);
	vector<librnary::LoopRegion> lrs;
	librnary::ExtractMultiLoopRegions(lrs, sst.RootSurface());
	EXPECT_EQ(lrs.size(), 1);
	EXPECT_EQ(librnary::ExtractBranches(lrs.front()), 4);
	EXPECT_EQ(librnary::ExtractUnpaired(lrs.front()), 6);
	EXPECT_EQ(librnary::ExtractSumAsymmetry(lrs.front()), 4);
	auto slengths = librnary::ExtractSegLengths(lrs.front());
	multiset<int> seglength(begin(slengths), end(slengths));
	EXPECT_EQ(seglength, multiset<int>({2,2,2,0}));
	vector<librnary::Surface> mlsurfs;
	librnary::ExtractMultiLoopSurfaces(mlsurfs, sst.RootSurface());
	EXPECT_EQ(mlsurfs.size(), 1);
	auto bszs = librnary::ExtractBranchSizes(mlsurfs.front());
	multiset<int> bsizes(begin(bszs), end(bszs));
	EXPECT_EQ(bsizes, multiset<int>({5,19,5,13}));
}