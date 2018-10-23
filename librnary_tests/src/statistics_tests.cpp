//
// Created by max on 5/30/16.
//

#include <gtest/gtest.h>

#include "secondary_structure.hpp"
#include "statistics.hpp"
#include "ss_enumeration.hpp"
#include "random.hpp"

using namespace std;

TEST(Statistics, SmallStem1) {
	auto ptrue = librnary::DotBracketToMatching("..(((...)))..");
	auto ppred = librnary::DotBracketToMatching("...((...))...");
	int tp = librnary::TruePositives(ptrue, ppred);
	int fp = librnary::FalsePositives(ptrue, ppred);
	int fn = librnary::FalseNegatives(ptrue, ppred);
	EXPECT_EQ(tp, 2);
	EXPECT_EQ(fp, 0);
	EXPECT_EQ(fn, 1);
}

TEST(Statistics, SmallStem2) {
	auto ptrue = librnary::DotBracketToMatching("(((((...)))))");
	auto ppred = librnary::DotBracketToMatching("...((...))...");
	int tp = librnary::TruePositives(ptrue, ppred);
	int fp = librnary::FalsePositives(ptrue, ppred);
	int fn = librnary::FalseNegatives(ptrue, ppred);
	EXPECT_EQ(tp, 2);
	EXPECT_EQ(fp, 0);
	EXPECT_EQ(fn, 3);
}

TEST(Statistics, InternalLoop1) {
	auto ptrue = librnary::DotBracketToMatching(".(((..((....)).)))");
	auto ppred = librnary::DotBracketToMatching(".(.(..((....)))..)");
	int tp = librnary::TruePositives(ptrue, ppred);
	int fp = librnary::FalsePositives(ptrue, ppred);
	int fn = librnary::FalseNegatives(ptrue, ppred);
	EXPECT_EQ(tp, 3);
	EXPECT_EQ(fp, 1);
	EXPECT_EQ(fn, 2);
	double sensitivity = 3.0 / 5.0, ppv = 3.0 / 4.0, f1score = 2.0 * ((ppv * sensitivity) / (ppv + sensitivity));
	EXPECT_DOUBLE_EQ(librnary::Sensitivity(ptrue, ppred), sensitivity);
	EXPECT_DOUBLE_EQ(librnary::PositivePredictiveValue(ptrue, ppred), ppv);
	EXPECT_DOUBLE_EQ(librnary::F1Score(ptrue, ppred), f1score);
}

TEST(Statistics, Bulge1) {
	auto ptrue = librnary::DotBracketToMatching("..(...((...))).");
	auto ppred = librnary::DotBracketToMatching("..(..(((...))))");
	int tp = librnary::TruePositives(ptrue, ppred);
	int fp = librnary::FalsePositives(ptrue, ppred);
	int fn = librnary::FalseNegatives(ptrue, ppred);
	EXPECT_EQ(tp, 2);
	EXPECT_EQ(fp, 2);
	EXPECT_EQ(fn, 1);
	double sensitivity = 2.0 / 3.0, ppv = 2.0 / 4.0, f1score = 2.0 * ((ppv * sensitivity) / (ppv + sensitivity));
	EXPECT_DOUBLE_EQ(librnary::Sensitivity(ptrue, ppred), sensitivity);
	EXPECT_DOUBLE_EQ(librnary::PositivePredictiveValue(ptrue, ppred), ppv);
	EXPECT_DOUBLE_EQ(librnary::F1Score(ptrue, ppred), f1score);
}

TEST(Statistics, FuzzTest) {
	const int CASES = 500, RNALEN = 20, TRIALS = 40;
	default_random_engine re = librnary::RandomEngineForTests();
	for (int tc = 0; tc < CASES; ++tc) {
		auto primary = librnary::RandomPrimary(re, RNALEN);
		unsigned trials = min((unsigned) librnary::BondPairs(primary).size(), (unsigned) TRIALS);
		auto ptrue = librnary::RandomMatching(primary, re, trials);
		auto ppred = librnary::RandomMatching(primary, re, trials);
		int correct_pairs = 0, true_pairs = 0, prediced_pairs = 0;
		for (int i = 0; i < RNALEN; ++i) {
			if (ptrue[i] != i) {
				++true_pairs;
			}
			if (ppred[i] != i) {
				++prediced_pairs;
			}
			if (ppred[i] != i && ptrue[i] == ppred[i]) {
				++correct_pairs;
			}
		}
		correct_pairs /= 2;
		true_pairs /= 2;
		prediced_pairs /= 2;
		int tp = librnary::TruePositives(ptrue, ppred);
		int fp = librnary::FalsePositives(ptrue, ppred);
		int fn = librnary::FalseNegatives(ptrue, ppred);
		EXPECT_EQ(tp, correct_pairs);
		EXPECT_EQ(tp + fn, true_pairs);
		EXPECT_EQ(tp + fp, prediced_pairs);
		double sensitivity = correct_pairs / double(true_pairs);
		double ppv = correct_pairs / double(prediced_pairs);
		double f1score = 2.0 * ((ppv * sensitivity) / (ppv + sensitivity));
		if (ppv + sensitivity <= 0.0) {
			f1score = 0.0;
		}
		EXPECT_DOUBLE_EQ(librnary::Sensitivity(ptrue, ppred), sensitivity);
		EXPECT_DOUBLE_EQ(librnary::PositivePredictiveValue(ptrue, ppred), ppv);
		EXPECT_DOUBLE_EQ(librnary::F1Score(ptrue, ppred), f1score);
	}
}

TEST(SlippageStatistics, SmallStem1) {
	auto ptrue = librnary::DotBracketToMatching("..(((...)))..");
	auto ppred = librnary::DotBracketToMatching("...((...))...");
	int tp = librnary::slippage::TruePositives(ptrue, ppred);
	int fp = librnary::slippage::FalsePositives(ptrue, ppred);
	int fn = librnary::slippage::FalseNegatives(ptrue, ppred);
	EXPECT_EQ(tp, 2);
	EXPECT_EQ(fp, 0);
	EXPECT_EQ(fn, 1);
}

TEST(SlippageStatistics, SmallStem2) {
	auto ptrue = librnary::DotBracketToMatching("(((((...)))))");
	auto ppred = librnary::DotBracketToMatching("...((...))...");
	int tp = librnary::slippage::TruePositives(ptrue, ppred);
	int fp = librnary::slippage::FalsePositives(ptrue, ppred);
	int fn = librnary::slippage::FalseNegatives(ptrue, ppred);
	EXPECT_EQ(tp, 2);
	EXPECT_EQ(fp, 0);
	EXPECT_EQ(fn, 3);
}

TEST(SlippageStatistics, InternalLoop1) {
	auto ptrue = librnary::DotBracketToMatching(".(((..((....)).)))");
	auto ppred = librnary::DotBracketToMatching(".(.(..((....)))..)");
	int tp = librnary::slippage::TruePositives(ptrue, ppred);
	int fp = librnary::slippage::FalsePositives(ptrue, ppred);
	int fn = librnary::slippage::FalseNegatives(ptrue, ppred);
	EXPECT_EQ(tp, 4);
	EXPECT_EQ(fp, 0);
	EXPECT_EQ(fn, 1);
	double sensitivity = tp / double(tp + fn), ppv = tp / double(tp + fp),
		f1score = 2.0 * ((ppv * sensitivity) / (ppv + sensitivity));
	EXPECT_DOUBLE_EQ(librnary::slippage::Sensitivity(ptrue, ppred), sensitivity);
	EXPECT_DOUBLE_EQ(librnary::slippage::PositivePredictiveValue(ptrue, ppred), ppv);
	EXPECT_DOUBLE_EQ(librnary::slippage::F1Score(ptrue, ppred), f1score);
}

TEST(SlippageStatistics, Bulge1) {
	auto ptrue = librnary::DotBracketToMatching("..(...((...))).");
	auto ppred = librnary::DotBracketToMatching("..(..(((...))))");
	int tp = librnary::slippage::TruePositives(ptrue, ppred);
	int fp = librnary::slippage::FalsePositives(ptrue, ppred);
	int fn = librnary::slippage::FalseNegatives(ptrue, ppred);
	EXPECT_EQ(tp, 3);
	EXPECT_EQ(fp, 1);
	EXPECT_EQ(fn, 0);
	double sensitivity = tp / double(tp + fn), ppv = tp / double(tp + fp),
		f1score = 2.0 * ((ppv * sensitivity) / (ppv + sensitivity));
	EXPECT_DOUBLE_EQ(librnary::slippage::Sensitivity(ptrue, ppred), sensitivity);
	EXPECT_DOUBLE_EQ(librnary::slippage::PositivePredictiveValue(ptrue, ppred), ppv);
	EXPECT_DOUBLE_EQ(librnary::slippage::F1Score(ptrue, ppred), f1score);
}

TEST(SlippageStatistics, FalseNeg) {
	auto ptrue = librnary::DotBracketToMatching(".(...)...");
	auto ppred = librnary::DotBracketToMatching("....(...)");
	int tp = librnary::slippage::TruePositives(ptrue, ppred);
	int fp = librnary::slippage::FalsePositives(ptrue, ppred);
	int fn = librnary::slippage::FalseNegatives(ptrue, ppred);
	EXPECT_EQ(tp, 0);
	EXPECT_EQ(fp, 1);
	EXPECT_EQ(fn, 1);
	double sensitivity = tp / double(tp + fn), ppv = tp / double(tp + fp),
		f1score = 2.0 * ((ppv * sensitivity) / (ppv + sensitivity));
	if (ppv + sensitivity == 0)
		f1score = 0;
	EXPECT_DOUBLE_EQ(librnary::slippage::Sensitivity(ptrue, ppred), sensitivity);
	EXPECT_DOUBLE_EQ(librnary::slippage::PositivePredictiveValue(ptrue, ppred), ppv);
	EXPECT_DOUBLE_EQ(librnary::slippage::F1Score(ptrue, ppred), f1score);
}

TEST(SlippageStatistics, FalseNeg2) {
	auto ptrue = librnary::DotBracketToMatching("..((...))");
	auto ppred = librnary::DotBracketToMatching("....(...)");
	int tp = librnary::slippage::TruePositives(ptrue, ppred);
	int fp = librnary::slippage::FalsePositives(ptrue, ppred);
	int fn = librnary::slippage::FalseNegatives(ptrue, ppred);
	EXPECT_EQ(tp, 0);
	EXPECT_EQ(fp, 1);
	EXPECT_EQ(fn, 2);
	double sensitivity = tp / double(tp + fn), ppv = tp / double(tp + fp),
		f1score = 2.0 * ((ppv * sensitivity) / (ppv + sensitivity));
	if (ppv + sensitivity == 0)
		f1score = 0;
	EXPECT_DOUBLE_EQ(librnary::slippage::Sensitivity(ptrue, ppred), sensitivity);
	EXPECT_DOUBLE_EQ(librnary::slippage::PositivePredictiveValue(ptrue, ppred), ppv);
	EXPECT_DOUBLE_EQ(librnary::slippage::F1Score(ptrue, ppred), f1score);
}

