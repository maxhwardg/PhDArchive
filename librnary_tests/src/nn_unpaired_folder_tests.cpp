//
// Created by max on 6/19/16.
//

#include <gtest/gtest.h>
#include "models/nn_unpaired_model.hpp"
#include <scorers/nn_scorer.hpp>
#include <ss_enumeration.hpp>
#include "folders/nn_unpaired_folder.hpp"
#include <folders/brute_folder.hpp>
#include "models/nn_affine_model.hpp"

#include "random.hpp"

using namespace std;

const string DATA_TABLE_PATH = "../../data_tables/";

TEST(NNUnpairedFolder, ZeroOnEmptyRNA) {
	librnary::NNUnpairedModel model(DATA_TABLE_PATH);
	librnary::NNUnpairedFolder folder(model);
	auto prim = librnary::StringToPrimary("");
	EXPECT_EQ(0, folder.Fold(prim));
}

TEST(NNUnpairedFolder, GGGAGCUCCAUUGCUGACCUGCAGGGAUCGAAUCCGUAUGGGGCACUUAU) {
	librnary::NNUnpairedModel model(DATA_TABLE_PATH);
	librnary::NNScorer<librnary::NNUnpairedModel> scorer(model);
	librnary::NNUnpairedFolder folder(model);
	auto prim = librnary::StringToPrimary("GGGAGCUCCAUUGCUGACCUGCAGGGAUCGAAUCCGUAUGGGGCACUUAU");
	int fold_mfe = folder.Fold(prim);
	auto fold_trace = folder.Traceback();
	librnary::SSTree sstree(fold_trace);
	scorer.SetRNA(prim);
	EXPECT_EQ(fold_mfe, scorer.ScoreExterior(sstree.RootSurface()));
}


TEST(NNUnpairedFolder, TRNALikeMultiWithLinearParams) {
	auto dt = librnary::LoadDatatable(DATA_TABLE_PATH);
	librnary::NNUnpairedModel model(DATA_TABLE_PATH);
	model.SetMLParams(dt->efn2a, dt->efn2c, dt->efn2b, 0, 999999);
	librnary::NNScorer<librnary::NNUnpairedModel> scorer(model);
	librnary::NNUnpairedFolder folder(model);
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

TEST(NNUnpairedFolder, InternalLoopWithRequiredLonelyPair) {
	librnary::NNUnpairedModel model(DATA_TABLE_PATH);
	librnary::NNScorer<librnary::NNUnpairedModel> scorer(model);
	librnary::NNUnpairedFolder folder(model);
	auto prim = librnary::StringToPrimary("CGUCAAUAACCGUAGGAGUAUGGGGUCUGUAAGCAAAAGUAUUGUAUCGGGUACGACAUUGGG");
	// Expected structure: .(((.....(((..(.(((((...((......))....))))).)..)))....)))......
	int fold_mfe = folder.Fold(prim);
	auto fold_trace = folder.Traceback();
	librnary::SSTree sstree(fold_trace);
	scorer.SetRNA(prim);
	EXPECT_EQ(fold_mfe, -97);
	EXPECT_EQ(fold_mfe, scorer.ScoreExterior(sstree.RootSurface()));
	folder.SetLonelyPairs(false);
	EXPECT_LT(fold_mfe, folder.Fold(prim));
	EXPECT_NE(fold_trace, folder.Traceback());
}


TEST(NNUnpairedFolder, NestedMulti_GUUCUGAGCAUCGCUCCUAGGCGCCGA) {
	librnary::NNUnpairedModel model(DATA_TABLE_PATH);
	model.SetMLParams(-92, -83, 124, 4.0, 19);
	librnary::NNScorer<librnary::NNUnpairedModel> scorer(model);
	librnary::NNUnpairedFolder folder(model);
	auto prim = librnary::StringToPrimary("GUUCUGAGCAUCGCUCCUAGGCGCCGA");
	// ..((...)((....)(...))(...))
	int fold_mfe = folder.Fold(prim);
	auto fold_trace = folder.Traceback();
	librnary::SSTree sstree(fold_trace);
	scorer.SetRNA(prim);
	EXPECT_EQ(fold_mfe, -628);
	EXPECT_EQ(fold_mfe, scorer.ScoreExterior(sstree.RootSurface()));
}


TEST(NNUnpairedFolder, NestedMultiWithUnpaired_AAAUCAAUUUCGGAUUGUUGCCGAGA) {
	librnary::NNUnpairedModel model(DATA_TABLE_PATH);
	model.SetMLParams(-350, -24, -51, 3.8, 2);
	librnary::NNScorer<librnary::NNUnpairedModel> scorer(model);
	librnary::NNUnpairedFolder folder(model);
	auto prim = librnary::StringToPrimary("AAAUCAAUUUCGGAUUGUUGCCGAGA");
	// ...(.(.(...).(...)).(...))
	int fold_mfe = folder.Fold(prim);
	auto fold_trace = folder.Traceback();
	librnary::SSTree sstree(fold_trace);
	scorer.SetRNA(prim);
	EXPECT_EQ(fold_mfe, -910);
	EXPECT_EQ(fold_mfe, scorer.ScoreExterior(sstree.RootSurface()));
}


TEST(NNUnpairedFolder, SensitiveToMultiUnpaired) {
	librnary::NNUnpairedModel model(DATA_TABLE_PATH);
	model.SetMLParams(-182, 110, -146, -1.2, 18);
	librnary::NNScorer<librnary::NNUnpairedModel> scorer(model);
	librnary::NNUnpairedFolder folder(model);
	auto prim = librnary::StringToPrimary("UAGCAGGUGCUGCCGAAACUUU");
	// .(......(...).(...)..)
	folder.SetMaxMulti(9);
	int fold_mfe = folder.Fold(prim);
	auto fold_trace = folder.Traceback();
	librnary::SSTree sstree(fold_trace);
	scorer.SetRNA(prim);
	EXPECT_EQ(fold_mfe, -1096);
	EXPECT_EQ(fold_mfe, scorer.ScoreExterior(sstree.RootSurface()));
	folder.SetMaxMulti(8);
	EXPECT_LE(fold_mfe, folder.Fold(prim));
}

// TODO: Tests for limited feature size.

//TEST(NNUnpairedFolder, BigMultiWithRequiredUnpaired) {
//	librnary::NNUnpairedModel model(DATA_TABLE_PATH);
//	librnary::NNScorer<librnary::NNUnpairedModel> scorer(model);
//	librnary::NNUnpairedFolder folder(model);
//	auto prim = librnary::StringToPrimary("AACCCAGUCAAGCGACCGACACUUGGCGCACUUCGUACAGUUCAGUGUGGACUGCGCCUUAGAGUUUUGCCC"
//											  "ACAUACAGCAGAGUUUUGCGGAUUGGGUAAUCUAUGGAUUAUAAUUCCUUCCUGGAGUGAAAUUUAC"
//											  "ACCCUUG");
//	// Expected structure:
//	// .((((((((..((((..(((.((.((((((.(((((((......))))))).))))))..)).)))))))........
//	// .((((....)))).))))))))......((((((...(((((.....)))))..))))))........
//	folder.SetMaxMulti(12);
//	int fold_mfe = folder.Fold(prim);
//	auto fold_trace = folder.Traceback();
//	librnary::SSTree sstree(fold_trace);
//	scorer.SetRNA(prim);
//	EXPECT_EQ(fold_mfe, -397);
//	EXPECT_EQ(fold_mfe, scorer.ScoreExterior(sstree.RootSurface()));
//	// At least 11 unpaired in a multi-loop are needed for MFE on this structure.
//	folder.SetMaxMulti(10);
//	EXPECT_LT(fold_mfe, folder.Fold(prim));
//	EXPECT_NE(fold_trace, folder.Traceback());
//}

