//
// Created by max on 5/29/16.
//

#include <gtest/gtest.h>
#include <stack>

#include "ss_tree.hpp"
#include "ss_enumeration.hpp"
#include "random.hpp"

using namespace std;

TEST(SSTree, TwoInExternal1) {
	auto pairing = librnary::DotBracketToMatching("(()..(..).)..(..)");
	librnary::SSTree rna_tree(pairing);
	auto root = rna_tree.RootId();
	EXPECT_EQ(rna_tree.NumChildren(root), 2);
	auto left_child = rna_tree.Child(root, 0);
	EXPECT_EQ(rna_tree.NumChildren(left_child), 2);
	auto right_child = rna_tree.Child(root, 1);
	EXPECT_EQ(rna_tree.NumChildren(right_child), 0);
	EXPECT_EQ(rna_tree.Parent(right_child), root);
	EXPECT_EQ(rna_tree.Parent(left_child), root);
	for (int i = 0; i < rna_tree.NumChildren(left_child); ++i) {
		EXPECT_EQ(rna_tree.Parent(rna_tree.Child(left_child, i)), left_child);
		EXPECT_EQ(rna_tree.NumChildren(rna_tree.Child(left_child, i)), 0);
	}

}

TEST(SSTree, Nested1) {
	auto pairing = librnary::DotBracketToMatching("...(..(..(..).())..)");
	librnary::SSTree rna_tree(pairing);
	auto root = rna_tree.RootId();
	EXPECT_EQ(rna_tree.NumChildren(root), 1);
	auto child = rna_tree.Child(root, 0);
	EXPECT_EQ(rna_tree.NumChildren(child), 1);
	auto childchild = rna_tree.Child(child, 0);
	EXPECT_EQ(rna_tree.NumChildren(childchild), 2);
}

TEST(SSTree, Nested2) {
	auto pairing = librnary::DotBracketToMatching("((..))");
	librnary::SSTree rna_tree(pairing);
	auto root = rna_tree.RootId();
	EXPECT_EQ(rna_tree.NumChildren(root), 1);
	auto child = rna_tree.Child(root, 0);
	EXPECT_EQ(rna_tree.Parent(child), root);
	EXPECT_EQ(rna_tree.NumChildren(child), 1);
	auto childchild = rna_tree.Child(child, 0);
	EXPECT_EQ(rna_tree.Parent(childchild), child);
	EXPECT_EQ(rna_tree.NumChildren(childchild), 0);
}

TEST(SSTree, FuzzDFSFindsAllPairings) {
	default_random_engine re = librnary::RandomEngineForTests();
	const int RNA_LEN = 50, CASES = 1000, TRIALS = 50;
	for (int tc = 0; tc < CASES; ++tc) {
		// Generate a random pairing.
		auto prim = librnary::RandomPrimary(re, RNA_LEN);
		auto pairing = librnary::RandomMatching(prim, re, TRIALS);
		int unpaired = 0;
		for (int i = 0; i < static_cast<int>(pairing.size()); ++i) {
			if (pairing[i] == i)
				++unpaired;
		}
		// DFS through the pairing.
		librnary::SSTree rna_tree(pairing);
		stack<librnary::SSTreeNodeId> s;
		s.push(rna_tree.RootId());
		// Check assumption about external loop.
		EXPECT_TRUE(rna_tree.IsExternalLoop(s.top()));
		int dfs_unpaired = 0;
		while (!s.empty()) {
			auto curr = s.top();
			s.pop();
			dfs_unpaired += rna_tree.Unpaired(curr);
			for (auto c : rna_tree.Children(curr)) {
				// Check assumption about external loop.
				EXPECT_FALSE(rna_tree.IsExternalLoop(c));
				// Mark pair locations as seen, and make sure they are valid.
				EXPECT_EQ(pairing[rna_tree.PairI(c)], rna_tree.PairJ(c));
				// Make sure the parent relationship is valid.
				EXPECT_EQ(rna_tree.Parent(c), curr);
				pairing[rna_tree.PairI(c)] = rna_tree.PairI(c);
				pairing[rna_tree.PairJ(c)] = rna_tree.PairJ(c);
				s.push(c);
			}
		}
		// Check all pair locations were seen.
		for (int i = 0; i < static_cast<int>(pairing.size()); ++i) {
			EXPECT_EQ(pairing[i], i);
		}
		// Check that the correct #unpaired was seen during dfs.
		EXPECT_EQ(unpaired, dfs_unpaired);
	}
}

TEST(Surface, FuzzDFSFindsAllPairings) {
	default_random_engine re = librnary::RandomEngineForTests();
	const int RNA_LEN = 50, CASES = 1000, TRIALS = 50;
	for (int tc = 0; tc < CASES; ++tc) {
		// Generate a random pairing.
		auto prim = librnary::RandomPrimary(re, RNA_LEN);
		auto pairing = librnary::RandomMatching(prim, re, TRIALS);
		int unpaired = 0;
		for (int i = 0; i < static_cast<int>(pairing.size()); ++i) {
			if (pairing[i] == i)
				++unpaired;
		}
		// DFS through the pairing.
		librnary::SSTree rna_tree(pairing);
		stack<librnary::Surface> s;
		s.push(rna_tree.RootSurface());
		// Check assumption about external loop.
		EXPECT_TRUE(s.top().IsExternalLoop());
		int dfs_unpaired = 0;
		while (!s.empty()) {
			auto curr = s.top();
			s.pop();
			dfs_unpaired += curr.Unpaired();
			for (auto &c : curr.Children()) {
				// Check assumption about external loop.
				EXPECT_FALSE(c.IsExternalLoop());
				// Mark pair locations as seen, and make sure they are valid.
				EXPECT_EQ(pairing[c.PairI()], c.PairJ());
				// Make sure the parent relationship is valid.
				EXPECT_EQ(c.Parent(), curr);
				pairing[c.PairI()] = c.PairI();
				pairing[c.PairJ()] = c.PairJ();
				s.push(c);
			}
		}
		// Check all pair locations were seen.
		for (int i = 0; i < static_cast<int>(pairing.size()); ++i) {
			EXPECT_EQ(pairing[i], i);
		}
		// Check that the correct #unpaired was seen during dfs.
		EXPECT_EQ(unpaired, dfs_unpaired);
	}
}