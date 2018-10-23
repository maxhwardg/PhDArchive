//
// Created by max on 8/30/17.
//

#ifndef RNARK_SCORER_UTILS_HPP
#define RNARK_SCORER_UTILS_HPP

#include "gtest/gtest.h"

#include "ss_tree.hpp"
#include "scorers/nn_scorer.hpp"

#include <stack>

template <typename Scorer>
void CheckTrace(Scorer scorer, librnary::SSTree tree) {
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
		EXPECT_EQ(stack_score, at.stacking);
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

#endif //RNARK_SCORER_UTILS_HPP
