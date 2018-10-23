//
// Created by max on 6/28/16.
//

#include <gtest/gtest.h>
#include "models/aalberts_model.hpp"
#include "scorers/aalberts_scorer.hpp"
#include "folders/aalberts_folder.hpp"
#include <folders/brute_folder.hpp>
#include "models/nn_unpaired_model.hpp"
#include "folders/nn_unpaired_folder.hpp"

#include "random.hpp"

using namespace std;

const string DATA_TABLE_PATH = "../../data_tables/";

TEST(AalbertsFolder, ZeroForEmptyRNA) {
	librnary::AalbertsModel model(DATA_TABLE_PATH);
	librnary::AalbertsFolder folder(model);
	auto prim = librnary::StringToPrimary("");
	int fold_mfe = folder.Fold(prim);
	auto fold_trace = folder.Traceback();
	EXPECT_EQ(fold_mfe, 0);
	EXPECT_EQ(fold_trace.size(), 0);
}

TEST(AalbertsFolder, CGCCCUGGCUCAUGAUGUCUCAGUGGCUCUACAUCGGGCU) {
	librnary::AalbertsModel model(DATA_TABLE_PATH);
	librnary::AalbertsScorer scorer(model);
	librnary::AalbertsFolder folder(model);
	folder.SetMaxALength(10);
	folder.SetMaxBLength(5);
	auto prim = librnary::StringToPrimary("CGCCCUGGCUCAUGAUGUCUCAGUGGCUCUACAUCGGGCU");
	int fold_mfe = folder.Fold(prim);
	auto fold_trace = folder.Traceback();
	auto ss_tree = librnary::SSTree(fold_trace);
	scorer.SetRNA(prim);
	EXPECT_EQ(fold_mfe, scorer.ScoreExterior(ss_tree.RootSurface()));
}

TEST(AalbertsFolder, CAUCACGAGCAGCUCCGCCUUUUGGACUCCUUUGUACUGGUGACGAUGAU) {
	librnary::AalbertsModel model(DATA_TABLE_PATH);
	librnary::AalbertsScorer scorer(model);
	librnary::AalbertsFolder folder(model);
	folder.SetMaxALength(10);
	folder.SetMaxBLength(5);
	auto prim = librnary::StringToPrimary("CAUCACGAGCAGCUCCGCCUUUUGGACUCCUUUGUACUGGUGACGAUGAU");
	int fold_mfe = folder.Fold(prim);
	auto fold_trace = folder.Traceback();
	auto ss_tree = librnary::SSTree(fold_trace);
	scorer.SetRNA(prim);
	EXPECT_EQ(fold_mfe, scorer.ScoreExterior(ss_tree.RootSurface()));
}

TEST(AalbertsFolder, GUCCCAGAGUGGGUGGUACCACGCGGGCUGACGAGCUUACGAAGCCUAGG) {
	librnary::AalbertsModel model(DATA_TABLE_PATH);
	librnary::AalbertsScorer scorer(model);
	librnary::AalbertsFolder folder(model);
	folder.SetMaxALength(10);
	folder.SetMaxBLength(5);
	auto prim = librnary::StringToPrimary("GUCCCAGAGUGGGUGGUACCACGCGGGCUGACGAGCUUACGAAGCCUAGG");
	int fold_mfe = folder.Fold(prim);
	auto fold_trace = folder.Traceback();
	auto ss_tree = librnary::SSTree(fold_trace);
	scorer.SetRNA(prim);
	EXPECT_EQ(fold_mfe, scorer.ScoreExterior(ss_tree.RootSurface()));
}

TEST(AalbertsFolder, GAGGCGGUAUGGAAGGCCGUUAACGACGGCUAAACCCCACGUAACAACAAGUGUCCAUUG) {
	librnary::AalbertsModel model(DATA_TABLE_PATH);
	librnary::AalbertsScorer scorer(model);
	librnary::AalbertsFolder folder(model);
	folder.SetMaxALength(10);
	folder.SetMaxBLength(5);
	auto prim = librnary::StringToPrimary("GAGGCGGUAUGGAAGGCCGUUAACGACGGCUAAACCCCACGUAACAACAAGUGUCCAUUG");
	int fold_mfe = folder.Fold(prim);
	auto fold_trace = folder.Traceback();
	auto ss_tree = librnary::SSTree(fold_trace);
	scorer.SetRNA(prim);
	EXPECT_EQ(fold_mfe, scorer.ScoreExterior(ss_tree.RootSurface()));
}

TEST(AalbertsFolder, GCGCCCGCGCGGACUUUGGCAUGGGUCGGCCAAGUUACGCCGUCGAAAGUAUGUCAGCAC) {
	librnary::AalbertsModel model(DATA_TABLE_PATH);
	librnary::AalbertsScorer scorer(model);
	librnary::AalbertsFolder folder(model);
	folder.SetMaxALength(10);
	folder.SetMaxBLength(5);
	auto prim = librnary::StringToPrimary("GCGCCCGCGCGGACUUUGGCAUGGGUCGGCCAAGUUACGCCGUCGAAAGUAUGUCAGCAC");
	int fold_mfe = folder.Fold(prim);
	auto fold_trace = folder.Traceback();
	auto ss_tree = librnary::SSTree(fold_trace);
	scorer.SetRNA(prim);
	EXPECT_EQ(fold_mfe, scorer.ScoreExterior(ss_tree.RootSurface()));
}

TEST(AalbertsFolder, AUUUCAGAACACUUCUGGCAGGGAAAUCCGGGUGAGGCCAUC) {
	librnary::AalbertsModel model(DATA_TABLE_PATH);
	model.SetMLParams(0, 46, 1, 1, 1);
	librnary::AalbertsScorer scorer(model);
	librnary::AalbertsFolder folder(model);
	folder.SetMaxALength(10);
	folder.SetMaxBLength(5);
	auto prim = librnary::StringToPrimary("AUUUCAGAACACUUCUGGCAGGGAAAUCCGGGUGAGGCCAUC");
	// Expecting ...((((((...))))))..(((...))).((((....))))
	int fold_mfe = folder.Fold(prim);
	auto fold_trace = folder.Traceback();
	auto ss_tree = librnary::SSTree(fold_trace);
	scorer.SetRNA(prim);
	EXPECT_EQ(fold_mfe, scorer.ScoreExterior(ss_tree.RootSurface()));
}


TEST(AalbertsFolder, GAAUAAGCGCCGUGGAAGGGGCUUGCUGAAUGGUCUGGAUGCCGCCCUAAGUAAGCCGCC) {
	librnary::AalbertsModel model(DATA_TABLE_PATH);
	librnary::AalbertsScorer scorer(model);
	librnary::AalbertsFolder folder(model);
	folder.SetMaxALength(10);
	folder.SetMaxBLength(5);
	auto prim = librnary::StringToPrimary("GAAUAAGCGCCGUGGAAGGGGCUUGCUGAAUGGUCUGGAUGCCGCCCUAAGUAAGCCGCC");
	int fold_mfe = folder.Fold(prim);
	auto fold_trace = folder.Traceback();
	auto ss_tree = librnary::SSTree(fold_trace);
	scorer.SetRNA(prim);
	model.SetRNA(prim);
	cout << model.MLInit(3, 1) + model.FlushCoax(9, 18, 19, 56) + model.ClosingThreeDangle(7, 57) << endl;
	EXPECT_EQ(fold_mfe, scorer.ScoreExterior(ss_tree.RootSurface()));
}


TEST(AalbertsFolder, MultiLoopOrder3FlushCoaxialStackInside) {
	librnary::AalbertsModel model(DATA_TABLE_PATH);
	librnary::AalbertsScorer scorer(model);
	librnary::AalbertsFolder folder(model);
	folder.SetMaxALength(10);
	folder.SetMaxBLength(5);
	auto prim = librnary::StringToPrimary("UGAUUUGAGCAAGUGUCUUGUCUAAUAUAACAC"
											  "CUCAGGGUCGUAGGAAACGCGACCCCCUUGGAGGCGCCCUUAGCGAAAGGCUCGCUAGCGUGUUGUA");
	int fold_mfe = folder.Fold(prim);
	auto fold_trace = folder.Traceback();
	auto ss_tree = librnary::SSTree(fold_trace);
	scorer.SetRNA(prim);
	EXPECT_EQ(fold_mfe, scorer.ScoreExterior(ss_tree.RootSurface()));
}

// TODO: Tests for limited feature size.
