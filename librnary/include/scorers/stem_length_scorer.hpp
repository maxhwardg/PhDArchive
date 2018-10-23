//
// Created by max on 3/19/18.
//

#ifndef RNARK_STEM_LENGTH_SCORER_HPP
#define RNARK_STEM_LENGTH_SCORER_HPP

#include "models/stem_length_model.hpp"
#include "scorers/nn_scorer.hpp"

namespace librnary {

/**
 * This class extends the classic NNModel scorer to work with the StemLengthModel.
 */
class StemLengthScorer: public NNScorer<StemLengthModel> {
protected:
	/// If true, single nucleotide bulges will not break a stack.
	bool allow_single_nt_bulges = false;
public:

	/// If true, single nucleotide bulges will not break a stack.
	bool AllowSingleNtBulges() const;

	/// If true, single nucleotide bulges will not break a stack.
	void SetAllowSingleNtBulges(bool to);

	/**
	 * Overrides the NNScorer ScoreInternal to include the stem length score.
	 * @param surf The super surface that encloses the structure to score.
	 * @return The free energy change.
	 */
	energy_t ScoreInternal(const Surface &surf) const override;

	SurfaceScore TraceInternal(const Surface &surf) const override;

	StemLengthScorer(StemLengthModel _em) : NNScorer(_em, true) {}

	StemLengthScorer(StemLengthModel _em, bool _stacking) : NNScorer(_em, _stacking) {}
};

/**
 * A convenince class that is like StemLengthScorer, but has AllowSingleNtBulges set to true by default.
 * Useful for templating reasons, such as in hotfuzz.
 */
class StemLengthBulgeScorer: public StemLengthScorer {
public:
	StemLengthBulgeScorer(StemLengthModel _em) : StemLengthScorer(_em, true) {
		allow_single_nt_bulges = true;
	}

	StemLengthBulgeScorer(StemLengthModel _em, bool _stacking) : StemLengthScorer(_em, _stacking) {
		allow_single_nt_bulges = true;
	}
};

}

#endif //RNARK_STEM_LENGTH_SCORER_HPP
