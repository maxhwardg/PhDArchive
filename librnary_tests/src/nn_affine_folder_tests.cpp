//
// Created by max on 8/6/17.
//

#include <gtest/gtest.h>

#include "folders/nn_affine_folder.hpp"
#include "scorers/nn_scorer.hpp"

using namespace std;

const string DATA_TABLE_PATH = "../../data_tables/";

TEST(NNAffineFolder, ZeroOnEmptyRNA) {
	librnary::NNAffineModel model(DATA_TABLE_PATH);
	librnary::NNAffineFolder folder(model);
	auto prim = librnary::StringToPrimary("");
	EXPECT_EQ(0, folder.Fold(prim));
}

TEST(NNAffineFolder, GGGAGCUCCAUUGCUGACCUGCAGGGAUCGAAUCCGUAUGGGGCACUUAU) {
	librnary::NNAffineModel model(DATA_TABLE_PATH);
	librnary::NNScorer<librnary::NNAffineModel> scorer(model);
	librnary::NNAffineFolder folder(model);
	auto prim = librnary::StringToPrimary("GGGAGCUCCAUUGCUGACCUGCAGGGAUCGAAUCCGUAUGGGGCACUUAU");
	int fold_mfe = folder.Fold(prim);
	auto fold_trace = folder.Traceback();
	librnary::SSTree sstree(fold_trace);
	scorer.SetRNA(prim);
	EXPECT_EQ(fold_mfe, scorer.ScoreExterior(sstree.RootSurface()));
}


TEST(NNAffineFolder, TRNALikeMultiWithLinearParams) {
	auto dt = librnary::LoadDatatable(DATA_TABLE_PATH);
	librnary::NNAffineModel model(DATA_TABLE_PATH);
	model.SetMLParams(dt->efn2a, dt->efn2c, dt->efn2b);
	librnary::NNScorer<librnary::NNAffineModel> scorer(model);
	librnary::NNAffineFolder folder(model);
	folder.SetLonelyPairs(false);
	auto prim = librnary::StringToPrimary("ACAAGAGAUAGCGAAGCUGCGCUGCCGGAUGGUUUUUCAUCCUAACUUUUGAGUUCUGUGCUUACCAUAGCGACC"
											  "CCACACGUC");
	// .((((((.(((((......)))))..((((((....))))))...)))))).((((((((.....))))).)))..........
	int fold_mfe = folder.Fold(prim);
	auto fold_trace = folder.Traceback();
	librnary::SSTree sstree(fold_trace);
	scorer.SetRNA(prim);
	EXPECT_EQ(fold_mfe, -204);
	EXPECT_EQ(fold_mfe, scorer.ScoreExterior(sstree.RootSurface()));
}

TEST(NNAffineFolder, MLInitBonus) {
	auto dt = librnary::LoadDatatable(DATA_TABLE_PATH);
	librnary::NNAffineModel model(DATA_TABLE_PATH);
	model.SetMLParams(-183, 11, 1);
	librnary::NNScorer<librnary::NNAffineModel> scorer(model);
	librnary::NNAffineFolder folder(model);
	auto prim = librnary::StringToPrimary("UUGGACGGCCUAACAGAACGUUUAUAGAAGGGUAUAUCCGCUAGAGC");
	// ..((.(.(.((...))(...))(....).)((....))).(....))
	int fold_mfe = folder.Fold(prim);
	auto fold_trace = folder.Traceback();
	librnary::SSTree sstree(fold_trace);
	scorer.SetRNA(prim);
	EXPECT_EQ(fold_mfe, -523);
	EXPECT_EQ(fold_mfe, scorer.ScoreExterior(sstree.RootSurface()));
}