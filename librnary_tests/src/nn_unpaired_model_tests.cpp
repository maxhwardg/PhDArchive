//
// Created by max on 24/08/16.
//

#include <gtest/gtest.h>
#include "models/nn_unpaired_model.hpp"

using namespace std;

const string DATA_TABLE_PATH = "../../data_tables/";

TEST(NNUnpairedModel, SimpleMulti) {
	librnary::NNUnpairedModel model(DATA_TABLE_PATH);
	auto sst = librnary::SSTree(librnary::DotBracketToMatching("(.(...)(...))"));
	model.SetMLParams(103, -45, 2, 1.2, 3);
	model.SetMLInitConstant(102);
	EXPECT_EQ(model.MLClosure(sst.RootSurface().Child(0)), 102 - 45 * 3 + 2);
}

TEST(NNUnpairedModel, SimpleMultiManyUnpaired) {
	librnary::NNUnpairedModel model(DATA_TABLE_PATH);
	auto sst = librnary::SSTree(librnary::DotBracketToMatching("(.(...)...(...))"));
	model.SetMLParams(102, -12, 2, 1.2, 2);
	model.SetMLBranchCost(-45);
	model.SetMLUnpairedPivot(3);
	EXPECT_EQ(model.MLClosure(sst.RootSurface().Child(0)),
			  102 - 45 * 3 + 3 * 2 + librnary::KCalToEnergy(1.2 * log(4.0 / 3)));
}

TEST(NNUnpairedModel, FourWayMulti) {
	librnary::NNUnpairedModel model(DATA_TABLE_PATH);
	auto sst = librnary::SSTree(librnary::DotBracketToMatching("((...).(...)...(...).)"));
	model.SetMLParams(99, -45, 2, 1.2, 3);
	model.SetMLLogMultiplier(2.1);
	EXPECT_EQ(model.MLClosure(sst.RootSurface().Child(0)),
			  99 - 45 * 4 + 3 * 2 + librnary::KCalToEnergy(2.1 * log(5.0 / 3)));
}


TEST(NNUnpairedModel, UnstrainedThreeWay) {
	librnary::StrainedUnpairedModel model(DATA_TABLE_PATH);
	auto sst = librnary::SSTree(librnary::DotBracketToMatching("(.(...)...(...).)"));
	model.SetMLParams(99, -45, 2, 1.2, 3);
	model.SetMLLogMultiplier(2.1);
	model.SetMLStrain(23);
	EXPECT_EQ(model.MLClosure(sst.RootSurface().Child(0)),
			  99 - 45 * 3 + 3 * 2 + librnary::KCalToEnergy(2.1 * log(5.0 / 3)));
}

TEST(NNUnpairedModel, StrainedThreeWay) {
	librnary::StrainedUnpairedModel model(DATA_TABLE_PATH);
	auto sst = librnary::SSTree(librnary::DotBracketToMatching("((...).(...))"));
	model.SetMLParams(99, -45, 2, 1.2, 3);
	model.SetMLLogMultiplier(2.1);
	model.SetMLStrain(23);
	EXPECT_EQ(model.MLClosure(sst.RootSurface().Child(0)),
			  99 - 45 * 3 + 2 + 23);
}