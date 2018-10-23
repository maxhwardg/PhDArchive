/**
 * Created by max on 3/19/18.
 * Contains tests related to the Stem Length Model, Scorer, and Folder.
 */

#include <gtest/gtest.h>
#include <folders/stem_length_folder.hpp>

#include "scorers/stem_length_scorer.hpp"
#include "paths.hpp"

namespace librnary {

TEST(StemLengthModel, Just8) {
	StemLengthModel em(DATA_TABLE_PATH);
	std::vector<energy_t> costs = {8};
	em.SetLengthCosts(costs);
	EXPECT_EQ(em.StemLengthCost(1), 8);
	EXPECT_EQ(em.StemLengthCost(2), 8);
	EXPECT_EQ(em.StemLengthCost(3), 8);
	EXPECT_EQ(em.StemLengthCosts(), costs);
}

TEST(StemLengthModel, Costs10and12) {
	StemLengthModel em(DATA_TABLE_PATH);
	std::vector<energy_t> costs = {10, 12};
	em.SetLengthCosts(costs);
	EXPECT_EQ(em.StemLengthCost(1), 10);
	EXPECT_EQ(em.StemLengthCost(2), 12);
	EXPECT_EQ(em.StemLengthCost(3), 12);
	EXPECT_EQ(em.StemLengthCost(100), 12);
	EXPECT_EQ(em.StemLengthCosts(), costs);
}


TEST(StemLengthScorer, LargerLeftBulgeWithLarge3PrimeHelix) {
	auto rna = librnary::StringToPrimary("GGAACCCCAAAGGGGUUU");
	//										    "GGAACCCCAAAGGGGUUU"
	auto match = librnary::DotBracketToMatching("(((.((((...)))))))");
	auto sst = librnary::SSTree(match);

	StemLengthModel em(DATA_TABLE_PATH, rna);
	std::vector<energy_t> costs = {10, 12, -6, -8};
	em.SetLengthCosts(costs);
	StemLengthScorer scorer(em);
	scorer.SetRNA(rna);

	NNAffineModel nnaffine_em(DATA_TABLE_PATH);
	NNScorer<NNAffineModel> nnscorer(nnaffine_em);
	nnscorer.SetRNA(rna);

	EXPECT_EQ(scorer.ScoreExterior(sst.RootSurface()), nnscorer.ScoreExterior(sst.RootSurface()) - 6 - 8);
}


TEST(StemLengthScorer, SimpleMulti) {
	auto rna = librnary::StringToPrimary("ACCGAACAUGGUUCCGAGCGUCU");
	//										    "ACCGAACAUGGUUCCGAGCGUCU"
	auto match = librnary::DotBracketToMatching("(..((((...)))).((...)))");
	auto sst = librnary::SSTree(match);

	StemLengthModel em(DATA_TABLE_PATH, rna);
	std::vector<energy_t> costs = {10, 12, -6};
	em.SetLengthCosts(costs);
	StemLengthScorer scorer(em);
	scorer.SetRNA(rna);

	NNAffineModel nnaffine_em(DATA_TABLE_PATH);
	NNScorer<NNAffineModel> nnscorer(nnaffine_em);
	nnscorer.SetRNA(rna);

	EXPECT_EQ(scorer.ScoreExterior(sst.RootSurface()), nnscorer.ScoreExterior(sst.RootSurface()) + 10 - 6 + 12);
}

TEST(StemLengthScorer, BulgesAndInternal) {
	auto rna = librnary::StringToPrimary("AGAGUGGGAACCUUAU");
	//										    "AGAGUGGGAACCUUAU"
	auto match = librnary::DotBracketToMatching(".(..(..(...)..))");
	auto sst = librnary::SSTree(match);

	StemLengthModel em(DATA_TABLE_PATH, rna);
	std::vector<energy_t> costs = {10, 12, -6};
	em.SetLengthCosts(costs);
	StemLengthScorer scorer(em);
	scorer.SetRNA(rna);

	NNAffineModel nnaffine_em(DATA_TABLE_PATH);
	NNScorer<NNAffineModel> nnscorer(nnaffine_em);
	nnscorer.SetRNA(rna);

	EXPECT_EQ(scorer.ScoreExterior(sst.RootSurface()), nnscorer.ScoreExterior(sst.RootSurface()) + 10 * 3);
}



TEST(StemLengthBulgeScorer, IgnoresSingleBulge1) {
	auto rna = librnary::StringToPrimary("UGGCGUUAC");
	auto match = librnary::DotBracketToMatching(".((...).)");
	auto sst = librnary::SSTree(match);

	StemLengthModel em(DATA_TABLE_PATH, rna);
	std::vector<energy_t> costs = {10, 12, -6};
	em.SetLengthCosts(costs);
	StemLengthScorer scorer(em);
	scorer.SetRNA(rna);
	scorer.SetAllowSingleNtBulges(true);


	NNAffineModel nnaffine_em(DATA_TABLE_PATH);
	NNScorer<NNAffineModel> nnscorer(nnaffine_em);
	nnscorer.SetRNA(rna);

	EXPECT_EQ(scorer.ScoreExterior(sst.RootSurface()), nnscorer.ScoreExterior(sst.RootSurface()) + 12);
}

TEST(StemLengthBulgeScorer, IgnoresSingleBulge2) {
	auto rna = librnary::StringToPrimary("AAGGGGUAAGCCCUCUUUCUAGGCAUUACGAGCGU");
	//											 AAGGGGUAAGCCCUCUUUCAGGCAUUACGAGCGU
	auto match = librnary::DotBracketToMatching(".(..((.(((...).))))..((..(...))).)");
	auto sst = librnary::SSTree(match);

	StemLengthModel em(DATA_TABLE_PATH, rna);
	std::vector<energy_t> costs = {10, 12, -6, 3, 6, 7};
	em.SetLengthCosts(costs);
	StemLengthScorer scorer(em);
	scorer.SetRNA(rna);
	scorer.SetAllowSingleNtBulges(true);


	NNAffineModel nnaffine_em(DATA_TABLE_PATH);
	NNScorer<NNAffineModel> nnscorer(nnaffine_em);
	nnscorer.SetRNA(rna);

	EXPECT_EQ(scorer.ScoreExterior(sst.RootSurface()), nnscorer.ScoreExterior(sst.RootSurface()) + 10 * 2 + 12 + 6);
}

TEST(StemLengthFolder, UGGCGUUC) {
	auto rna = librnary::StringToPrimary("UGGCGUUC");
	auto match = librnary::DotBracketToMatching(".((...))");

	StemLengthModel em(DATA_TABLE_PATH, rna);
	std::vector<energy_t> costs = {2, -61, 69, 18};
	em.SetLengthCosts(costs);
	StemLengthFolder folder(em);

	EXPECT_EQ(folder.Fold(rna), -17);
	EXPECT_EQ(folder.Traceback(), match);
}

}