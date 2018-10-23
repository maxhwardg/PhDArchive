//
// Created by max on 8/11/16.
//

#include <gtest/gtest.h>
#include "models/average_asym_model.hpp"

using namespace std;

const string DATA_TABLE_PATH = "../../data_tables/";

TEST(AverageAsymmetryModel, SimpleMultiLoop) {
	auto sst = librnary::SSTree(librnary::DotBracketToMatching("(..(...).(...))"));
	librnary::AverageAsymmetryModel model(DATA_TABLE_PATH);
	model.SetMLParams(10, -3, 2, 0, 1.5, 0.3);
	model.SetMLInit(11);
	EXPECT_EQ(model.MLClosure(sst.RootSurface().Child(0)),
			  11 + -3 * 3 + 2 * 3 + static_cast<int>(round(10 * (4.0 / 3) * 0.3)));
}

TEST(AverageAsymmetryModel, LargeAsymmetryMulti) {
	auto sst = librnary::SSTree(librnary::DotBracketToMatching("(.....(...).(...).)"));
	librnary::AverageAsymmetryModel model(DATA_TABLE_PATH);
	model.SetMLParams(10, -3, 2, 0, 1.5, 0.3);
	model.SetMLBranchCost(-4);
	model.SetMLMaxAvgAsymmetry(1.9);
	EXPECT_EQ(model.MLClosure(sst.RootSurface().Child(0)),
			  10 + -4 * 3 + 2 * 7 + static_cast<int>(round(10 * 1.9 * 0.3)));
}

TEST(AverageAsymmetryModel, LargeAsymmetryMultiNoLimitOnAsym) {
	auto sst = librnary::SSTree(librnary::DotBracketToMatching("(.....(...).(...).)"));
	librnary::AverageAsymmetryModel model(DATA_TABLE_PATH);
	model.SetMLParams(10, -3, 2, 0, 1.5, 0.3);
	model.SetMLBranchCost(-4);
	model.SetMLMaxAvgAsymmetry(3.0);
	EXPECT_EQ(model.MLClosure(sst.RootSurface().Child(0)),
			  10 + -4 * 3 + 2 * 7 + static_cast<int>(round(10 * (8.0 / 3) * 0.3)));
}

TEST(AverageAsymmetryModel, 5WayBranchMultiLoop) {
	auto sst = librnary::SSTree(librnary::DotBracketToMatching("((...)....(...).(...).((...))..)"));
	librnary::AverageAsymmetryModel model(DATA_TABLE_PATH);
	model.SetMLParams(10, -3, 2, 0, 1.5, 0.3);
	model.SetMLAsymmetryCoeff(1.2);
	model.SetMLUnpairedCost(-5);
	EXPECT_EQ(model.MLClosure(sst.RootSurface().Child(0)),
			  10 + -3 * 5 - 5 * 8 + static_cast<int>(round(10 * 1.5 * 1.2)));
}

TEST(AverageAsymmetryModel, 5WayBranchMultiLoopNoMaxAsym) {
	auto sst = librnary::SSTree(librnary::DotBracketToMatching("((...)....(...).(...).((...))..)"));
	librnary::AverageAsymmetryModel model(DATA_TABLE_PATH);
	model.SetMLParams(10, -3, 2, 0, 1.5, 0.3);
	model.SetMLAsymmetryCoeff(1.2);
	model.SetMLUnpairedCost(-5);
	model.SetMLMaxAvgAsymmetry(10.0);
	EXPECT_EQ(model.MLClosure(sst.RootSurface().Child(0)),
			  10 + -3 * 5 - 5 * 8 + static_cast<int>(round(10 * 10.0 / 5 * 1.2)));
}