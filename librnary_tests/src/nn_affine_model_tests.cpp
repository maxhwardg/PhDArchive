//
// Created by max on 7/31/17.
//

#include <gtest/gtest.h>

#include "models/nn_affine_model.hpp"

using namespace std;

const string DATA_TABLE_PATH = "../../data_tables/";

TEST(NNAffineModel, SimpleMulti) {
	librnary::NNAffineModel model(DATA_TABLE_PATH);
	auto sst = librnary::SSTree(librnary::DotBracketToMatching("(.(...)(...))"));
	model.SetMLParams(103, -45, 2);
	model.SetMLInitCost(102);
	EXPECT_EQ(model.MLClosure(sst.RootSurface().Child(0)), 102 - 45 * 3 + 2);
}

TEST(NNAffineModel, SimpleMultiManyUnpaired) {
	librnary::NNAffineModel model(DATA_TABLE_PATH);
	auto sst = librnary::SSTree(librnary::DotBracketToMatching("(.(...)...(...))"));
	model.SetMLParams(102, -12, 2);
	model.SetMLBranchCost(-45);
	EXPECT_EQ(model.MLClosure(sst.RootSurface().Child(0)), 102 - 45 * 3 + 4 * 2);
}

TEST(NNAffineModel, FourWayMulti) {
	librnary::NNAffineModel model(DATA_TABLE_PATH);
	auto sst = librnary::SSTree(librnary::DotBracketToMatching("((...).(...)...(...).)"));
	model.SetMLParams(99, -45, 10);
	model.SetMLUnpairedCost(2);
	EXPECT_EQ(model.MLClosure(sst.RootSurface().Child(0)), 99 - 45 * 4 + 5 * 2);
}