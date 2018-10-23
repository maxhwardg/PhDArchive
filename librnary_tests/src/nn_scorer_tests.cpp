//
// Created by max on 6/14/16.
//

#include <gtest/gtest.h>
#include <scorers/nn_scorer.hpp>
#include <models/nn_model.hpp>
#include "models/nn_affine_model.hpp"
#include <ss_enumeration.hpp>

#include "random.hpp"
#include "scorer_utils.hpp"

using namespace std;

const string DATA_TABLE_PATH = "../../data_tables/";



TEST(NNScorer, OneNucleotide) {
	librnary::NNAffineModel model(DATA_TABLE_PATH);
	librnary::NNScorer<librnary::NNAffineModel> scorer(model);
	auto prim = librnary::StringToPrimary("U");
	auto match = librnary::DotBracketToMatching(".");
	scorer.SetRNA(prim);
	librnary::SSTree tree(match);
	auto struc = librnary::LoadStructure(prim, match);
	auto dt = librnary::LoadDatatable(DATA_TABLE_PATH);
	EXPECT_EQ(scorer.ScoreExterior(tree.RootSurface()), 0);
	EXPECT_EQ(scorer.ScoreExterior(tree.RootSurface()), librnary::RunEFN2WithSimpleMulti(*dt, *struc));
	CheckTrace(scorer, tree);
}

TEST(NNScorer, SmallFavourableStem) {
	librnary::NNAffineModel model(DATA_TABLE_PATH);
	librnary::NNScorer<librnary::NNAffineModel> scorer(model);
	auto prim = librnary::StringToPrimary("ACGGAUAGGUCC");
	auto match = librnary::DotBracketToMatching("..(((....)))");
	scorer.SetRNA(prim);
	librnary::SSTree tree(match);
	auto struc = librnary::LoadStructure(prim, match);
	auto dt = librnary::LoadDatatable(DATA_TABLE_PATH);
	EXPECT_EQ(scorer.ScoreExterior(tree.RootSurface()), -5);
	EXPECT_EQ(scorer.ScoreExterior(tree.RootSurface()), librnary::RunEFN2WithSimpleMulti(*dt, *struc));
	CheckTrace(scorer, tree);
}

TEST(NNScorer, SmallRecursiveMultiLoop) {
	librnary::NNAffineModel model(DATA_TABLE_PATH);
	librnary::NNScorer<librnary::NNAffineModel> scorer(model);
	auto prim = librnary::StringToPrimary("GUACGUUGACUUUGUUCCGAGUUUUAUCCAGAAUUCCUUCAU");
	auto match = librnary::DotBracketToMatching("..((((((...).(...)).)(...).)..)(.....))...");
	scorer.SetRNA(prim);
	librnary::SSTree tree(match);
	auto struc = librnary::LoadStructure(prim, match);
	auto dt = librnary::LoadDatatable(DATA_TABLE_PATH);
	EXPECT_EQ(scorer.ScoreExterior(tree.RootSurface()), 470);
	EXPECT_EQ(scorer.ScoreExterior(tree.RootSurface()), librnary::RunEFN2WithSimpleMulti(*dt, *struc));
	CheckTrace(scorer, tree);
}

TEST(NNScorer, GUGUFExternalLooplushCoax) {
	// This test is a breaking case for RNAstructure 5.8.
	// The bug has been manually fixed in the source included with rnark.
	librnary::NNAffineModel model(DATA_TABLE_PATH);
	librnary::NNScorer<librnary::NNAffineModel> scorer(model);
	auto prim = librnary::StringToPrimary("UCUGAGUAAAUUGCUACGCG");
	auto match = librnary::DotBracketToMatching("(....)((...).......)");
	scorer.SetRNA(prim);
	librnary::SSTree tree(match);
	auto struc = librnary::LoadStructure(prim, match);
	auto dt = librnary::LoadDatatable(DATA_TABLE_PATH);
	EXPECT_EQ(scorer.ScoreExterior(tree.RootSurface()), 179);
	EXPECT_EQ(scorer.ScoreExterior(tree.RootSurface()), librnary::RunEFN2WithSimpleMulti(*dt, *struc));
	CheckTrace(scorer, tree);
}

TEST(NNScorer, UGUGMultiLoopCoax) {
	// This test is a breaking case for RNAstructure 5.8.
	// The bug has been manually fixed in the source included with rnark.
	librnary::NNAffineModel model(DATA_TABLE_PATH);
	librnary::NNScorer<librnary::NNAffineModel> scorer(model);
	auto prim = librnary::StringToPrimary("GCGCAGGUCAAGUAGCAUGU");
	auto match = librnary::DotBracketToMatching("......((...)(.....))");
	scorer.SetRNA(prim);
	librnary::SSTree tree(match);
	auto struc = librnary::LoadStructure(prim, match);
	auto dt = librnary::LoadDatatable(DATA_TABLE_PATH);
	EXPECT_EQ(scorer.ScoreExterior(tree.RootSurface()), 204);
	EXPECT_EQ(scorer.ScoreExterior(tree.RootSurface()), librnary::RunEFN2WithSimpleMulti(*dt, *struc));
	CheckTrace(scorer, tree);
}

TEST(NNScorer, LargeComplexRandom1) {
	librnary::NNAffineModel model(DATA_TABLE_PATH);
	librnary::NNScorer<librnary::NNAffineModel> scorer(model);
	auto prim = librnary::StringToPrimary("GCCAGCCGUGAUACCCCUUUGGCUACCUUAGUGUGCGCCCAUACCGCUCGGAGUUUUGCUCAUCAACAAUAGAAUA"
											  "UUAAGCGGCAAGCGCCGCGGUCGAAGUGAGCUGGAUGUGACCACUGGGCGGUCACCGGCUUCCCGGGAGUGA"
											  "AGCAUGCAGUGUUUAUGCCGGAUCACUCCGGGCUCGGUAAUGUUUACAGUCGGAUAGCAAAUUUCCAGCGG"
											  "AACAGGAGGCACCUUAUCAAUGCCGGCUCCACCGUGUCCCUGAUAGGGGUCUCUCUAUCUUAAAAGACGGCA"
											  "CCCCAACAAAUGGUCCGCGUCUCUAGCAGUUUGAGUCUGGUCGAACAUCCGAAAACACUAUCGCAGGUCUGU"
											  "CCUACCGAUGGCUUAAGGCCCAUCAAUGGCCUACGGUUGGUCGUAGUAGCUUAAUCCACACACAGAGACUC"
											  "CGCGAUGGUGCCGAGAAGACUAGG");
	auto match = librnary::DotBracketToMatching("((.....))..(.(((....)..(((.(((((((....(.((.((.....(.(((((.((..(((..."
													".(...)).((..(...).)).).((....))))..))...)).((((((..(.....))))).("
													"(....)(...).))))..).)((.(.(.....(...(....))))(...).)....)(.(..."
													".).)))).((.(...)(..(((..........)....))....)((..(...))..(.((..."
													".).)..))((...).)).).)(...)))))((.....(.(....).((.(....)...)(.(."
													"(...)..)((.((...((...)(.....).))))...)..((...).)))))((....).)))."
													")..).(...)).)..)..(.((...).).(...(..(..(..(...))..).)(...))..).)"
													"..)....))");
	scorer.SetRNA(prim);
	librnary::SSTree tree(match);
	auto struc = librnary::LoadStructure(prim, match);
	auto dt = librnary::LoadDatatable(DATA_TABLE_PATH);
	EXPECT_EQ(scorer.ScoreExterior(tree.RootSurface()), 3503);
	EXPECT_EQ(scorer.ScoreExterior(tree.RootSurface()), librnary::RunEFN2WithSimpleMulti(*dt, *struc));
	CheckTrace(scorer, tree);
}

TEST(NNScorer, LargeComplexRandom2) {
	librnary::NNAffineModel model(DATA_TABLE_PATH);
	librnary::NNScorer<librnary::NNAffineModel> scorer(model);
	auto prim = librnary::StringToPrimary("UCUGGAAGUCCACGACAAGUCUGAUAAGCACAUUCACUCUAAUGCGUGAACGCCGCCCCAACAGGUUGAGGUA"
											  "CUCUUUCUUAGACCGGGUUGUUCUAAGUCCAUUGGUGAGAGGACUGGUGAAUUGUAAACUACGGGUUC"
											  "UACACACGAAUAAUAGAUAGUCAACCUCUCAUUGUGAAAGGGAAUGGUCGAACGAUUCGUGUAGAUGG"
											  "AAUGCUAUCCCACCAUCCAAACGCCGAAGACUUGGUUACGUGCUUGUCGAAUAAACCCCCUAGCAUAU"
											  "CGAUAUAACCUAACGCGCUUAUCAUGAUUGCAUACUGUGAAGGGGGGGAGCCGUCGUAUGCCAGCAUG"
											  "UAGGGGUGUCGCCGUCCCUUUUUCAUAGAUACCGUUUCAUGUGUGCCGCCUUAGUCCGCAGGUUUCGU"
											  "GGACAUAUGUAUAACGGGUCAUUGCUUUUAAUACGUUUGAUCUUCGGCCUCACAUCAACGCAUAUUUU"
											  "GCGUUCAAAUACGUAGACUGAGAGAUCGUG");
	auto match = librnary::DotBracketToMatching(".(..((..(.....)((..(((.((..(.(...)...)...))).).(...(......)..))."
													".).((.(((...)((...)(.(...((((..(...)(...((....).).((.(..(.."
													".((.(...(.((((....)...)(.(.(...(.....)))(...).....))))).).)..("
													"((......)).....(....)))(..((...)......))(..((..(....)..).((..."
													"..((.((.......).).).....)).).).))..)....)((.(....).)..).).))..("
													"(...)).)..(((.(((...(.(....)..)).)((.((((((....))...))((...)."
													"((..(...))(..(((.(((((...).)...).....)((...)..)(....).))).)))."
													"))))))).(...)(....).).).))...(...)).)..)).).)(.((((.(...))..))"
													"(...))..)).)))..");
	scorer.SetRNA(prim);
	librnary::SSTree tree(match);
	auto struc = librnary::LoadStructure(prim, match);
	auto dt = librnary::LoadDatatable(DATA_TABLE_PATH);
	EXPECT_EQ(scorer.ScoreExterior(tree.RootSurface()), 4165);
	EXPECT_EQ(scorer.ScoreExterior(tree.RootSurface()), librnary::RunEFN2WithSimpleMulti(*dt, *struc));
	CheckTrace(scorer, tree);
}

TEST(NNScorer, LargeMFE) {
	librnary::NNAffineModel model(DATA_TABLE_PATH);
	librnary::NNScorer<librnary::NNAffineModel> scorer(model);
	auto prim = librnary::StringToPrimary("UGCCGUCAAUAACCGUAGGAGUAUGGGGUCUGUAAGCAAAAGUAUUGUAUCGGGUACGACAUUGGGGAACAUGGA"
											  "GCUCCUUUCCACGUUGAAAGAAGUAGUCACAAUCCUUCUGGGCAAAGGGUAAGAGAAGACUGCAGGUAU"
											  "GCAUCUUUGAAUCCUAAGAGACUUUGGCCCGUCGUGUUAAUUUAACGGCAGUUAU");
	auto match = librnary::DotBracketToMatching("((((((.((((((((..((.((..(((((((....((((.....)))).(((....)))....(("
													"(((...(((.(((.(((((......))))).)))..)))...)))))........(((((.."
													"..(((((.((((....)))))))))..)))))...))))))).))))..)).))))..)).))"
													")))).....");
	scorer.SetRNA(prim);
	librnary::SSTree tree(match);
	auto struc = librnary::LoadStructure(prim, match);
	auto dt = librnary::LoadDatatable(DATA_TABLE_PATH);
	EXPECT_EQ(scorer.ScoreExterior(tree.RootSurface()), -491);
	EXPECT_EQ(scorer.ScoreExterior(tree.RootSurface()), librnary::RunEFN2WithSimpleMulti(*dt, *struc));
	CheckTrace(scorer, tree);
}
