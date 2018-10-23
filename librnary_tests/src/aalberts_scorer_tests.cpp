//
// Created by max on 6/27/16.
//

#include <gtest/gtest.h>
#include "scorers/aalberts_scorer.hpp"
#include "models/nn_unpaired_model.hpp"
#include <ss_enumeration.hpp>

#include "random.hpp"

using namespace std;

const string DATA_TABLE_PATH = "../../data_tables/";

void CheckTrace(librnary::AalbertsScorer scorer, librnary::SSTree tree) {
	auto em = scorer.GetEnergyModel();
	std::stack<librnary::SurfaceScore> stk;
	stk.push(scorer.TraceExterior(tree.RootSurface()));
	while (!stk.empty()) {
		auto at = stk.top();
		stk.pop();
		librnary::energy_t stack_score = 0;
		for (const auto &stack_ft : at.stacks) {
			stack_score += stack_ft.Score(em, at.surface.PairI(), at.surface.PairJ());
		}
		// Different logic for external vs internal multi-loops.
		// TODO: Fix this different logic when OptimalMLConfig is fixed to split closure/stacking scores.
		if (at.surface.IsExternalLoop()) {
			EXPECT_EQ(stack_score, at.stacking);
		} else if (at.subsurfaces.size() > 1) {
			// Add closure cost for multi-loops.
			stack_score += em.MLInit(at.ml_closure_features["N"], at.ml_closure_features["M"]);
			EXPECT_EQ(stack_score, at.ml_closure);
		} else { // If it is a one or two loop, no stacking score expected.
			EXPECT_EQ(stack_score, 0);
		}
		if (at.surface.IsExternalLoop()) {
			EXPECT_EQ(scorer.ScoreExterior(at.surface), at.recursive_score);
		} else {
			EXPECT_EQ(scorer.ScoreInternal(at.surface), at.recursive_score);
		}
		for (const auto &child : at.subsurfaces) {
			stk.push(*child);
		}
	}
}

TEST(AalbertsScorer, BasicTest1) {
	librnary::AalbertsModel model(DATA_TABLE_PATH);
	librnary::AalbertsScorer scorer(model);
	auto prim = librnary::StringToPrimary("GAGAAACAGAAACC");
	auto match = librnary::DotBracketToMatching("(.(...).(...))");
	auto tree = librnary::SSTree(match);
	auto rsurf = tree.RootSurface();
	auto mlsurf = rsurf.Child(0);
	model.SetRNA(prim);
	scorer.SetRNA(prim);
	EXPECT_EQ(get<1>(scorer.OptimalMLConfig(mlsurf)) + get<0>(scorer.OptimalMLConfig(mlsurf)),
			  model.MLInit(5, 1) + model.ThreeDangle(2, 6) + model.FlushCoax(0, 13, 8, 12));
	CheckTrace(scorer, tree);
}

TEST(AalbertsScorer, BasicTest2) {
	librnary::AalbertsModel model(DATA_TABLE_PATH);
	librnary::AalbertsScorer scorer(model);
	auto prim = librnary::StringToPrimary("GAGAAACAGAAACAC");
	auto match = librnary::DotBracketToMatching("(.(...).(...).)");
	auto tree = librnary::SSTree(match);
	auto rsurf = tree.RootSurface();
	auto mlsurf = rsurf.Child(0);
	model.SetRNA(prim);
	scorer.SetRNA(prim);
	EXPECT_EQ(get<1>(scorer.OptimalMLConfig(mlsurf)) + get<0>(scorer.OptimalMLConfig(mlsurf)),
			  model.MLInit(4, 1) + model.ThreeDangle(8, 12) + model.MismatchCoax(2, 6, 0, 14));
	CheckTrace(scorer, tree);
}

TEST(AalbertsScorer, BasicTest3) {
	librnary::AalbertsModel model(DATA_TABLE_PATH);
	librnary::AalbertsScorer scorer(model);
	auto prim = librnary::StringToPrimary("GAGAAACAUAAAAAU");
	auto match = librnary::DotBracketToMatching("(.(...).(...).)");
	auto tree = librnary::SSTree(match);
	auto rsurf = tree.RootSurface();
	auto mlsurf = rsurf.Child(0);
	model.SetRNA(prim);
	scorer.SetRNA(prim);
	EXPECT_EQ(get<1>(scorer.OptimalMLConfig(mlsurf)) + get<0>(scorer.OptimalMLConfig(mlsurf)),
			  model.MLInit(4, 1) + model.ThreeDangle(8, 12) + model.MismatchCoax(2, 6, 0, 14));
	CheckTrace(scorer, tree);
}

TEST(AalbertsScorer, BasicTest4) {
	librnary::AalbertsModel model(DATA_TABLE_PATH);
	librnary::AalbertsScorer scorer(model);
	auto prim = librnary::StringToPrimary("GACGAAACACUAAAAGAU");
	auto match = librnary::DotBracketToMatching("(..(...)..(...)..)");
	auto tree = librnary::SSTree(match);
	auto rsurf = tree.RootSurface();
	auto mlsurf = rsurf.Child(0);
	model.SetRNA(prim);
	scorer.SetRNA(prim);
	EXPECT_EQ(get<1>(scorer.OptimalMLConfig(mlsurf)) + get<0>(scorer.OptimalMLConfig(mlsurf)),
			  model.MLInit(9, 3) + model.ThreeDangle(3, 7) + model.Mismatch(10, 14) + model.ClosingThreeDangle(0, 17));
	CheckTrace(scorer, tree);
}

TEST(AalbertsScorer, GUCCCAGAGUGGGUGGUACCACGCGGGCUGACGAGCUUACGAAGCCUAGG) {
	librnary::AalbertsModel model(DATA_TABLE_PATH);
	librnary::AalbertsScorer scorer(model);
	auto prim = librnary::StringToPrimary("GUCCCAGAGUGGGUGGUACCACGCGGGCUGACGAGCUUACGAAGCCUAGG");
	auto match = librnary::DotBracketToMatching("...((...((((......))))..(((((..((......)).))))).))");
	auto tree = librnary::SSTree(match);
	auto rsurf = tree.RootSurface();
	auto mlsurf = rsurf.Child(0).Child(0);
	model.SetRNA(prim);
	scorer.SetRNA(prim);
	EXPECT_EQ(get<1>(scorer.OptimalMLConfig(mlsurf)) + get<0>(scorer.OptimalMLConfig(mlsurf)),
			  model.ThreeDangle(8, 21) + model.MismatchCoax(4, 48, 24, 46) + model.MLInit(7, 1));
	CheckTrace(scorer, tree);
}

TEST(AalbertsScorer, LessThanOrEqNNScorerStacking) {
	// If we take the optimal stacking output from NNScorer and add the corresponding MLInit from Aalberts,
	// then it should always be >= the Aalberts scorer stacking score.
	// Mostly they will be equal, but sometimes +ve coaxial stacks can be favourable in an Aalberts model.
	librnary::AalbertsModel amodel(DATA_TABLE_PATH);
	librnary::NNUnpairedModel umodel(DATA_TABLE_PATH);
	librnary::AalbertsScorer ascorer(amodel);
	librnary::NNScorer<librnary::NNUnpairedModel> nnscorer(umodel);
	auto re = librnary::RandomEngineForTests();
	const int CASES = 1000, RNA_SIZE = 50, TRIALS = 100;
	for (int tc = 0; tc < CASES; ++tc) {
		auto prim = librnary::RandomPrimary(re, RNA_SIZE);
		auto match = librnary::RandomMatching(prim, re, TRIALS);
		amodel.SetRNA(prim);
		ascorer.SetRNA(prim);
		nnscorer.SetRNA(prim);
		librnary::SSTree sstree(match);
		CheckTrace(ascorer, sstree);
		stack<librnary::Surface> s;
		s.push(sstree.RootSurface());
		while (!s.empty()) {
			auto surf = s.top();
			s.pop();
			// Check multi-loop.
			if (!surf.IsExternalLoop() && surf.NumChildren() > 1) {
				auto ml_conf_trace = nnscorer.TraceMLConfig(surf);
				// +1 for the closing branch.
				int a = surf.Unpaired() + surf.NumChildren() + 1, b = surf.NumChildren() + 1;
				for (const auto &stack : std::get<0>(ml_conf_trace)) {
					if (stack.t == librnary::StackType::MMCX) {
						a -= 2; // Lose 4 gain 2.
						b -= 2;
					} else if (stack.t == librnary::StackType::FLUSHCX) {
						// a -= 0; Lose 2 gain 2.
						b -= 2;
					}
				}
				EXPECT_LE(get<0>(ascorer.OptimalMLConfig(surf)) + get<1>(ascorer.OptimalMLConfig(surf)),
						  amodel.MLInit(a, b) + std::get<0>(nnscorer.OptimalMLConfig(surf)));
			}
			// Add sub-surfaces.
			for (auto &ss : surf.Children()) {
				s.push(ss);
			}
		}
	}
}

// Tests if the same as NNScorer when using a equivalent model.
TEST(AalbertsScorer, SameAsNNScorer) {
	librnary::AalbertsModel amodel(DATA_TABLE_PATH);
	librnary::NNUnpairedModel umodel(DATA_TABLE_PATH);
	// Force models to be the same.
	amodel.SetMLParams(0, 0, 1, 1, 1);
	umodel.SetMLParams(0, 0, 0, 0, 99999);
	librnary::AalbertsScorer ascorer(amodel);
	librnary::NNScorer<librnary::NNUnpairedModel> nnscorer(umodel);
	auto re = librnary::RandomEngineForTests();
	const int CASES = 1000, RNA_SIZE = 50, TRIALS = 100;
	for (int tc = 0; tc < CASES; ++tc) {
		auto prim = librnary::RandomPrimary(re, RNA_SIZE);
		auto match = librnary::RandomMatching(prim, re, TRIALS);
		amodel.SetRNA(prim);
		ascorer.SetRNA(prim);
		nnscorer.SetRNA(prim);
		librnary::SSTree sstree(match);
		CheckTrace(ascorer, sstree);
		stack<librnary::Surface> s;
		s.push(sstree.RootSurface());
		while (!s.empty()) {
			auto surf = s.top();
			s.pop();
			// Check multi-loop.
			if (!surf.IsExternalLoop() && surf.NumChildren() > 1) {
				auto nn_conf = nnscorer.OptimalMLConfig(surf);
				auto aal_conf = ascorer.OptimalMLConfig(surf);
				EXPECT_EQ(get<0>(aal_conf) + get<1>(aal_conf), get<0>(nn_conf) + get<1>(nn_conf));
			}
			// Add sub-surfaces.
			for (auto &ss : surf.Children()) {
				s.push(ss);
			}
		}
	}
}