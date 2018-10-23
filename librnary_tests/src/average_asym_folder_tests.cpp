//
// Created by max on 8/11/16.
//

#include <gtest/gtest.h>
#include "folders/average_asym_folder.hpp"
#include "folders/nn_unpaired_folder.hpp"
#include "scorers/average_asym_scorer.hpp"

#include "random.hpp"

using namespace std;

const string DATA_TABLE_PATH = "../../data_tables/";

TEST(AverageAsymmetryFolder, ACAUGACACUAAAGGCGC) {
	librnary::AverageAsymmetryModel model(DATA_TABLE_PATH);
	model.SetMLParams(148, -96, 150, 0, 2.02028, -1.84906);
	librnary::AverageAsymmetryFolder folder(model);
	// The bug required these settings.
	folder.SetMaxMLBranches(3);
	folder.SetMaxMLNonClosingAsym(8);
	folder.SetUnpairedGap(8);
	auto prim = librnary::StringToPrimary("ACAUGACACUAAAGGCGC");
	librnary::SSTree sst(librnary::DotBracketToMatching(".(.((...)(...)).)."));
	librnary::AverageAsymmetryScorer scorer(model);
	scorer.SetRNA(prim);
	int mfe = folder.Fold(prim);
	EXPECT_EQ(mfe, scorer.ScoreExterior(sst.RootSurface()));
	librnary::SSTree mfesst(folder.Traceback());
	EXPECT_EQ(scorer.ScoreExterior(mfesst.RootSurface()), scorer.ScoreExterior(sst.RootSurface()));
}

TEST(AverageAsymmetryFolder, CGAGUGCCUAGUGC) {
	librnary::AverageAsymmetryModel model(DATA_TABLE_PATH);
	model.SetMLParams(-173, -71, -160, 0, 4.23311, 0.53386);
	librnary::AverageAsymmetryFolder folder(model);
	folder.SetMaxMLBranches(4);
	folder.SetMaxMLNonClosingAsym(6);
	folder.SetUnpairedGap(8);
	auto prim = librnary::StringToPrimary("CGAGUGCCUAGUGC");
	librnary::SSTree sst(librnary::DotBracketToMatching(".(.(...)(...))"));
	librnary::AverageAsymmetryScorer scorer(model);
	int mfe = folder.Fold(prim);
	scorer.SetRNA(prim);
	EXPECT_EQ(mfe, scorer.ScoreExterior(sst.RootSurface()));
	librnary::SSTree mfesst(folder.Traceback());
	EXPECT_EQ(scorer.ScoreExterior(mfesst.RootSurface()), scorer.ScoreExterior(sst.RootSurface()));
}

TEST(AverageAsymmetryFolder, CGUAGUGGACCCGAU) {
	librnary::AverageAsymmetryModel model(DATA_TABLE_PATH);
	model.SetMLParams(325, -190, -138, 0, 8.02713, 1.4503);
	librnary::AverageAsymmetryFolder folder(model);
	folder.SetMaxMLBranches(3);
	folder.SetMaxMLNonClosingAsym(7);
	folder.SetUnpairedGap(7);
	auto prim = librnary::StringToPrimary("CGUAGUGGACCCGAU");
	auto true_match = librnary::DotBracketToMatching(".((...)(...)..)");
	librnary::SSTree sst(true_match);

	librnary::AverageAsymmetryScorer scorer(model);
	int mfe = folder.Fold(prim);
	scorer.SetRNA(prim);
	EXPECT_EQ(mfe, scorer.ScoreExterior(sst.RootSurface()));
	librnary::SSTree mfesst(folder.Traceback());
	EXPECT_EQ(scorer.ScoreExterior(mfesst.RootSurface()), scorer.ScoreExterior(sst.RootSurface()));

}

TEST(AverageAsymmetryFolder, CUUCUGGGCAAAGGGUA) {
	librnary::AverageAsymmetryModel model(DATA_TABLE_PATH);
	model.SetMLParams(-280, 111, -161, 0, 0.393389, -0.0345322);
	librnary::AverageAsymmetryFolder folder(model);
	folder.SetMaxMLBranches(3);
	folder.SetMaxMLNonClosingAsym(7);
	folder.SetUnpairedGap(7);
	auto prim = librnary::StringToPrimary("CUUCUGGGCAAAGGGUA");
	auto true_match = librnary::DotBracketToMatching(".(.(...)(...)...)");
	librnary::SSTree sst(true_match);

	librnary::AverageAsymmetryScorer scorer(model);
	int mfe = folder.Fold(prim);
	scorer.SetRNA(prim);
	EXPECT_EQ(mfe, scorer.ScoreExterior(sst.RootSurface()));
	librnary::SSTree mfesst(folder.Traceback());
	EXPECT_EQ(scorer.ScoreExterior(mfesst.RootSurface()), scorer.ScoreExterior(sst.RootSurface()));

}

TEST(AverageAsymmetryFolder, UCCCCCCCUGGCUAUCUUGA) {
	librnary::AverageAsymmetryModel model(DATA_TABLE_PATH);
	model.SetMLParams(462, -24, -83, 0, 8.87734, -0.824442);
	librnary::AverageAsymmetryFolder folder(model);
	folder.SetMaxMLBranches(3);
	folder.SetMaxMLNonClosingAsym(7);
	folder.SetUnpairedGap(7);
	auto prim = librnary::StringToPrimary("UCCCCCCCUGGCUAUCUUGA");
	auto true_match = librnary::DotBracketToMatching("(....(...)...(...).)");
	librnary::SSTree sst(true_match);

	librnary::AverageAsymmetryScorer scorer(model);
	int mfe = folder.Fold(prim);
	scorer.SetRNA(prim);
	EXPECT_EQ(mfe, scorer.ScoreExterior(sst.RootSurface()));
	librnary::SSTree mfesst(folder.Traceback());
	clog << librnary::MatchingToDotBracket(folder.Traceback()) << endl;
	EXPECT_EQ(scorer.ScoreExterior(mfesst.RootSurface()), scorer.ScoreExterior(sst.RootSurface()));

}

TEST(AverageAsymmetryFolder, DoNotCrashOnEmptyRNA) {
	librnary::AverageAsymmetryModel model(DATA_TABLE_PATH);
	librnary::AverageAsymmetryFolder folder(model);
	auto prim = librnary::StringToPrimary("");
	EXPECT_NO_FATAL_FAILURE(folder.Fold(prim));
}

TEST(AverageAsymmetryFolder, DoNotCrashOnUUGA) {
	librnary::AverageAsymmetryModel model(DATA_TABLE_PATH);
	librnary::AverageAsymmetryFolder folder(model);
	auto prim = librnary::StringToPrimary("UUGA");
	EXPECT_NO_FATAL_FAILURE(folder.Fold(prim));
}


TEST(AverageAsymmetryFolder, AGACCGCAGAUCCAGAUGC) {
	librnary::AverageAsymmetryModel model(DATA_TABLE_PATH);
	librnary::NNUnpairedModel umodel(DATA_TABLE_PATH);
	int a = 500, b = -200, c = -20;
	model.SetMLParams(a, b, c, 0, 1.0, 0.0);
	umodel.SetMLParams(a, b, c, 0, 99999);
	librnary::AverageAsymmetryFolder folder(model);
	folder.SetMaxMLBranches(4);
	folder.SetMaxMLNonClosingAsym(15);
	folder.SetUnpairedGap(10);
	librnary::NNUnpairedFolder ufolder(umodel);
	auto prim = librnary::StringToPrimary("AGACCGCAGAUCCAGAUGC");
	int mfe = folder.Fold(prim);
	int umfe = ufolder.Fold(prim);
	clog << librnary::PrimaryToString(prim) << endl;
	EXPECT_EQ(mfe, umfe);
	clog << librnary::MatchingToDotBracket(ufolder.Traceback()) << endl;
}

TEST(AverageAsymmetryFolder, GGAAACGAAACGAAACGAAACC) {
	librnary::AverageAsymmetryModel model(DATA_TABLE_PATH);
	librnary::NNUnpairedModel umodel(DATA_TABLE_PATH);
	int a = 500, b = -200, c = -20;
	model.SetMLParams(a, b, c, 0, 1.0, 0.0);
	umodel.SetMLParams(a, b, c, 0, 99999);
	librnary::AverageAsymmetryFolder folder(model);
	folder.SetMaxMLBranches(5);
	folder.SetMaxMLNonClosingAsym(5);
	folder.SetUnpairedGap(5);
	librnary::NNUnpairedFolder ufolder(umodel);
	// Expected fold = ((...)(...)(...)(...))
	auto prim = librnary::StringToPrimary("GGAAACGAAACGAAACGAAACC");
	int mfe = folder.Fold(prim);
	int umfe = ufolder.Fold(prim);
	clog << librnary::PrimaryToString(prim) << endl;
	EXPECT_EQ(mfe, umfe);
	clog << librnary::MatchingToDotBracket(ufolder.Traceback()) << endl;
	clog << librnary::MatchingToDotBracket(folder.Traceback()) << endl;
	librnary::NNScorer<librnary::NNUnpairedModel> scorer(umodel);
	scorer.SetRNA(prim);
	clog << scorer.TraceExterior(librnary::SSTree(ufolder.Traceback()).RootSurface()).Describe(' ', 0) << endl;
}

TEST(AverageAsymmetryFolder, GGAAACGAAACGAAACAAAAAAC) {
	librnary::AverageAsymmetryModel model(DATA_TABLE_PATH);
	librnary::NNUnpairedModel umodel(DATA_TABLE_PATH);
	int a = 500, b = -200, c = -20;
	model.SetMLParams(a, b, c, 0, 1.0, 0.0);
	umodel.SetMLParams(a, b, c, 0, 99999);
	librnary::AverageAsymmetryFolder folder(model);
	folder.SetMaxMLBranches(4);
	folder.SetMaxMLNonClosingAsym(12);
	folder.SetUnpairedGap(6);
	librnary::NNUnpairedFolder ufolder(umodel);
	// Expected fold = ((...)(...)(...)......)
	auto prim = librnary::StringToPrimary("GGAAACGAAACGAAACAAAAAAC");
	int mfe = folder.Fold(prim);
	int umfe = ufolder.Fold(prim);

	EXPECT_EQ(mfe, umfe);

	librnary::NNScorer<librnary::NNUnpairedModel> scorer(umodel);
	scorer.SetRNA(prim);
	EXPECT_EQ(scorer.ScoreExterior(librnary::SSTree(folder.Traceback()).RootSurface()),
			  scorer.ScoreExterior(librnary::SSTree(ufolder.Traceback()).RootSurface()));

	clog << librnary::MatchingToDotBracket(ufolder.Traceback()) << endl;
	clog << librnary::MatchingToDotBracket(folder.Traceback()) << endl;
}

TEST(AverageAsymmetryFolder, GGAAACGAAACAAAAAAGAAACGAAACC) {
	librnary::AverageAsymmetryModel model(DATA_TABLE_PATH);
	librnary::NNUnpairedModel umodel(DATA_TABLE_PATH);
	int a = 500, b = -200, c = -20;
	model.SetMLParams(a, b, c, 0, 1.0, 0.0);
	umodel.SetMLParams(a, b, c, 0, 99999);
	librnary::AverageAsymmetryFolder folder(model);
	folder.SetMaxMLBranches(5);
	folder.SetMaxMLNonClosingAsym(12);
	folder.SetUnpairedGap(7);
	librnary::NNUnpairedFolder ufolder(umodel);
	// Expected fold = ((...)(...)......(...)(...)))
	auto prim = librnary::StringToPrimary("GGAAACGAAACAAAAAAGAAACGAAACC");
	int mfe = folder.Fold(prim);
	int umfe = ufolder.Fold(prim);

	EXPECT_EQ(mfe, umfe);


	librnary::NNScorer<librnary::NNUnpairedModel> scorer(umodel);
	scorer.SetRNA(prim);
	EXPECT_EQ(scorer.ScoreExterior(librnary::SSTree(folder.Traceback()).RootSurface()),
			  scorer.ScoreExterior(librnary::SSTree(ufolder.Traceback()).RootSurface()));

	clog << librnary::MatchingToDotBracket(ufolder.Traceback()) << endl;
	clog << librnary::MatchingToDotBracket(folder.Traceback()) << endl;
}

TEST(AverageAsymmetryFolder, GGAAACGAAACAAAAAAGAAACGAAACC_UnderLimitBranches) {
	librnary::AverageAsymmetryModel model(DATA_TABLE_PATH);
	librnary::NNUnpairedModel umodel(DATA_TABLE_PATH);
	int a = 500, b = -200, c = -20;
	model.SetMLParams(a, b, c, 0, 1.0, 0.0);
	umodel.SetMLParams(a, b, c, 0, 99999);
	librnary::AverageAsymmetryFolder folder(model);
	folder.SetMaxMLBranches(4);
	folder.SetMaxMLNonClosingAsym(12);
	folder.SetUnpairedGap(7);
	librnary::NNUnpairedFolder ufolder(umodel);
	// Expected fold = ((...)(...)......(...)(...)))
	auto prim = librnary::StringToPrimary("GGAAACGAAACAAAAAAGAAACGAAACC");
	int mfe = folder.Fold(prim);
	int umfe = ufolder.Fold(prim);

	EXPECT_GT(mfe, umfe);


	librnary::NNScorer<librnary::NNUnpairedModel> scorer(umodel);
	scorer.SetRNA(prim);
	EXPECT_GT(scorer.ScoreExterior(librnary::SSTree(folder.Traceback()).RootSurface()),
			  scorer.ScoreExterior(librnary::SSTree(ufolder.Traceback()).RootSurface()));

	clog << librnary::MatchingToDotBracket(ufolder.Traceback()) << endl;
	clog << librnary::MatchingToDotBracket(folder.Traceback()) << endl;
}

TEST(AverageAsymmetryFolder, GGAAACGAAACAAAAAAGAAACGAAACC_UnderLimitGap) {
	librnary::AverageAsymmetryModel model(DATA_TABLE_PATH);
	librnary::NNUnpairedModel umodel(DATA_TABLE_PATH);
	int a = 500, b = -200, c = -20;
	clog << a << " " << b << " " << c << endl;
	model.SetMLParams(a, b, c, 0, 1.0, 0.0);
	umodel.SetMLParams(a, b, c, 0, 99999);
	librnary::AverageAsymmetryFolder folder(model);
	folder.SetMaxMLBranches(5);
	folder.SetMaxMLNonClosingAsym(12);
	folder.SetUnpairedGap(5);
	librnary::NNUnpairedFolder ufolder(umodel);
	// Expected fold = ((...)(...)......(...)(...)))
	auto prim = librnary::StringToPrimary("GGAAACGAAACAAAAAAGAAACGAAACC");
	int mfe = folder.Fold(prim);
	int umfe = ufolder.Fold(prim);

	EXPECT_GT(mfe, umfe);


	librnary::NNScorer<librnary::NNUnpairedModel> scorer(umodel);
	scorer.SetRNA(prim);
	EXPECT_GT(scorer.ScoreExterior(librnary::SSTree(folder.Traceback()).RootSurface()),
			  scorer.ScoreExterior(librnary::SSTree(ufolder.Traceback()).RootSurface()));

	clog << librnary::MatchingToDotBracket(ufolder.Traceback()) << endl;
	clog << librnary::MatchingToDotBracket(folder.Traceback()) << endl;
}

TEST(AverageAsymmetryFolder, GGAAACGAAACAAAAAAGAAACGAAACC_UnderLimitAsym) {
	librnary::AverageAsymmetryModel model(DATA_TABLE_PATH);
	librnary::NNUnpairedModel umodel(DATA_TABLE_PATH);
	int a = 500, b = -200, c = -20;
	model.SetMLParams(a, b, c, 0, 1.0, 0.0);
	umodel.SetMLParams(a, b, c, 0, 99999);
	librnary::AverageAsymmetryFolder folder(model);
	folder.SetMaxMLBranches(5);
	folder.SetMaxMLNonClosingAsym(11);
	folder.SetUnpairedGap(6);
	librnary::NNUnpairedFolder ufolder(umodel);
	// Expected fold = ((...)(...)......(...)(...)))
	auto prim = librnary::StringToPrimary("GGAAACGAAACAAAAAAGAAACGAAACC");
	int mfe = folder.Fold(prim);
	int umfe = ufolder.Fold(prim);

	EXPECT_GT(mfe, umfe);


	librnary::NNScorer<librnary::NNUnpairedModel> scorer(umodel);
	scorer.SetRNA(prim);
	EXPECT_GT(scorer.ScoreExterior(librnary::SSTree(folder.Traceback()).RootSurface()),
			  scorer.ScoreExterior(librnary::SSTree(ufolder.Traceback()).RootSurface()));

	clog << librnary::MatchingToDotBracket(ufolder.Traceback()) << endl;
	clog << librnary::MatchingToDotBracket(folder.Traceback()) << endl;
}

TEST(AverageAsymmetryFolder, GAAAGAAACGAAACGAAACAAAC) {
	librnary::AverageAsymmetryModel model(DATA_TABLE_PATH);
	librnary::NNUnpairedModel umodel(DATA_TABLE_PATH);
	int a = 500, b = -200, c = -20;
	model.SetMLParams(a, b, c, 0, 1.0, 0.0);
	umodel.SetMLParams(a, b, c, 0, 99999);
	librnary::AverageAsymmetryFolder folder(model);
	folder.SetMaxMLBranches(5);
	folder.SetMaxMLNonClosingAsym(6);
	folder.SetUnpairedGap(4);
	librnary::NNUnpairedFolder ufolder(umodel);
	// Expected fold = (...(...)(...)(...)...)
	auto prim = librnary::StringToPrimary("GAAAGAAACGAAACGAAACAAAC");
	int mfe = folder.Fold(prim);
	int umfe = ufolder.Fold(prim);

	EXPECT_EQ(mfe, umfe);


	librnary::NNScorer<librnary::NNUnpairedModel> scorer(umodel);
	scorer.SetRNA(prim);
	EXPECT_EQ(scorer.ScoreExterior(librnary::SSTree(folder.Traceback()).RootSurface()),
			  scorer.ScoreExterior(librnary::SSTree(ufolder.Traceback()).RootSurface()));

	clog << librnary::MatchingToDotBracket(ufolder.Traceback()) << endl;
	clog << librnary::MatchingToDotBracket(folder.Traceback()) << endl;
}

TEST(AverageAsymmetryFolder, GAAAGAAACGAAACGAAACAAACTooSmallSumAsym) {
	librnary::AverageAsymmetryModel model(DATA_TABLE_PATH);
	librnary::NNUnpairedModel umodel(DATA_TABLE_PATH);
	int a = 500, b = -200, c = -20;
	model.SetMLParams(a, b, c, 0, 1.0, 0.0);
	umodel.SetMLParams(a, b, c, 0, 99999);
	librnary::AverageAsymmetryFolder folder(model);
	folder.SetMaxMLBranches(5);
	folder.SetMaxMLNonClosingAsym(5);
	folder.SetUnpairedGap(4);
	librnary::NNUnpairedFolder ufolder(umodel);
	// Expected fold = (...(...)(...)(...)...)
	auto prim = librnary::StringToPrimary("GAAAGAAACGAAACGAAACAAAC");
	int mfe = folder.Fold(prim);
	int umfe = ufolder.Fold(prim);

	EXPECT_GT(mfe, umfe);


	librnary::NNScorer<librnary::NNUnpairedModel> scorer(umodel);
	scorer.SetRNA(prim);
	EXPECT_GT(scorer.ScoreExterior(librnary::SSTree(folder.Traceback()).RootSurface()),
			  scorer.ScoreExterior(librnary::SSTree(ufolder.Traceback()).RootSurface()));

	clog << librnary::MatchingToDotBracket(ufolder.Traceback()) << endl;
	clog << librnary::MatchingToDotBracket(folder.Traceback()) << endl;
}

TEST(AverageAsymmetryFolder, MaxNonCLosingAsymmetryRight) {
	librnary::AverageAsymmetryModel model(DATA_TABLE_PATH);
	librnary::NNUnpairedModel umodel(DATA_TABLE_PATH);
	int a = 500, b = -200, c = -20;
	model.SetMLParams(a, b, c, 0, 1.0, 0.0);
	umodel.SetMLParams(a, b, c, 0, 99999);
	librnary::AverageAsymmetryFolder folder(model);
	folder.SetMaxMLBranches(4);
	folder.SetUnpairedGap(6);
	folder.SetStacking(false);
	librnary::NNUnpairedFolder ufolder(umodel);
	ufolder.SetStacking(false);
	// Expected fold = ((...)......(...)(...))
	auto prim = librnary::StringToPrimary("GGAAACAAAAAAGAAACGAAACC");
	int mfe = folder.Fold(prim);
	int umfe = ufolder.Fold(prim);

	EXPECT_EQ(mfe, umfe);

	librnary::NNScorer<librnary::NNUnpairedModel> scorer(umodel);
	scorer.SetRNA(prim);
	EXPECT_EQ(scorer.ScoreExterior(librnary::SSTree(folder.Traceback()).RootSurface()),
			  scorer.ScoreExterior(librnary::SSTree(ufolder.Traceback()).RootSurface()));

	clog << librnary::MatchingToDotBracket(folder.Traceback()) << endl;
	clog << librnary::MatchingToDotBracket(ufolder.Traceback()) << endl;

	folder.SetMaxMLNonClosingAsym(11);

	// Check it fails when too few non closing asymmetry is allowed.

	EXPECT_GT(folder.Fold(prim), umfe);

	clog << librnary::MatchingToDotBracket(folder.Traceback()) << endl;
}

TEST(AverageAsymmetryFolder, MaxNonCLosingAsymmetryLeft) {
	librnary::AverageAsymmetryModel model(DATA_TABLE_PATH);
	librnary::NNUnpairedModel umodel(DATA_TABLE_PATH);
	int a = 500, b = -200, c = -20;
	model.SetMLParams(a, b, c, 0, 1.0, 0.0);
	umodel.SetMLParams(a, b, c, 0, 99999);
	librnary::AverageAsymmetryFolder folder(model);
	folder.SetMaxMLBranches(4);
	folder.SetUnpairedGap(6);
	folder.SetStacking(false);
	librnary::NNUnpairedFolder ufolder(umodel);
	ufolder.SetStacking(false);
	// Expected fold = ((...)(...)......(...))
	auto prim = librnary::StringToPrimary("GGAAACGAAACAAAAAAGAAACC");
	int mfe = folder.Fold(prim);
	int umfe = ufolder.Fold(prim);

	EXPECT_EQ(mfe, umfe);

	librnary::NNScorer<librnary::NNUnpairedModel> scorer(umodel);
	scorer.SetRNA(prim);
	EXPECT_EQ(scorer.ScoreExterior(librnary::SSTree(folder.Traceback()).RootSurface()),
			  scorer.ScoreExterior(librnary::SSTree(ufolder.Traceback()).RootSurface()));

	clog << librnary::MatchingToDotBracket(folder.Traceback()) << endl;
	clog << librnary::MatchingToDotBracket(ufolder.Traceback()) << endl;

	folder.SetMaxMLNonClosingAsym(11);

	// Check it fails when too few non closing asymmetry is allowed.

	EXPECT_GT(folder.Fold(prim), umfe);

	clog << librnary::MatchingToDotBracket(folder.Traceback()) << endl;
}