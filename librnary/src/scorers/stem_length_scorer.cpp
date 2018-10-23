//
// Created by max on 3/19/18.
//

#include "scorers/stem_length_scorer.hpp"

namespace {

class NotImplemented: public std::logic_error {
public:
	NotImplemented() : std::logic_error("Function not yet implemented") {};
};

}

namespace librnary {

bool StemLengthScorer::AllowSingleNtBulges() const {
	return allow_single_nt_bulges;
}

void StemLengthScorer::SetAllowSingleNtBulges(bool to) {
	allow_single_nt_bulges = to;
}

energy_t StemLengthScorer::ScoreInternal(const librnary::Surface &surf) const {
	assert(surf.PairI() >= 0 && static_cast<int>(em.RNA().size()) > surf.PairJ());
	int length = 1;
	energy_t stacking_energy = 0;
	auto curr = surf;
	while (curr.NumChildren() == 1 && curr.Unpaired() <= 1) {
		if (curr.Unpaired() == 1 && !allow_single_nt_bulges) {
			break;
		}
		auto child = curr.Child(0);
		stacking_energy += em.TwoLoop(curr.PairI(), child.PairI(), child.PairJ(), curr.PairJ());
		curr = child;
		++length;
	}
	return em.StemLengthCost(length) + NNScorer::ScoreInternal(curr) + stacking_energy;
}

SurfaceScore StemLengthScorer::TraceInternal(const librnary::Surface &surf) const {
	throw NotImplemented();
}
}

