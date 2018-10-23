//
// Created by max on 8/2/17.
//

#ifndef RNARK_IBF_MULTILOOP_RMSEFSCORE_HPP
#define RNARK_IBF_MULTILOOP_RMSEFSCORE_HPP


#include <vector>
#include <iostream>
#include <set>

#include "read_cts.hpp"
#include "energy.hpp"
#include "ss_tree.hpp"
#include "multi_loop.hpp"
#include "statistics.hpp"
#include "parallel.hpp"
#include "vector_types.hpp"
#include "IBF_multiloop.hpp"

namespace librnary {

/**
 * A modification of IBFMultiLoop that uses the RMSE of set average fscores.
 * Takes the RMSE of average f-score errors for each group. Average f-score error is 1-fscore.
 * @tparam ParamSetT Parameter set type.
 * @tparam ModelT Energy model type.
 * @tparam ScorerT Energy model socrer type.
 */
template<typename ParamSetT, typename ModelT, typename ScorerT>
class IBFMultiLoopRMSEFscore : public IBFMultiLoop<ParamSetT, ModelT, ScorerT> {

protected:

	long FindBestParams(const std::vector<ParamSetT> &params, const ModelT &model) override {
		librnary::parallel_transform(params, this->param_scores, [=](const ParamSetT &pset) {
			auto local_model = model;
			pset.LoadInto(local_model);

			double sum_sq_error = 0;

			for (size_t ctg = 0; ctg < this->cts.size(); ++ctg) {
				VE true_energies = this->true_multi_base_energies[ctg];
				auto min_base_energies = true_energies;

				V<int> min_energy_choice(this->cts[ctg].size(), -1);

				for (size_t i = 0; i < this->cts[ctg].size(); ++i) {
					for (const auto &mi : this->true_multi_info[ctg][i]) {
						true_energies[i] += mi.MLClosure(local_model);
					}
				}

				std::vector<energy_t> min_energies = true_energies;
				for (size_t i = 0; i < this->cts[ctg].size(); ++i) {
					for (size_t j = 0; j < this->false_multi_info[ctg][i].size(); ++j) {
						energy_t e = this->false_multi_base_energies[ctg][i][j];
						min_base_energies[i] = std::min(min_base_energies[i], e);
						for (const auto &mi : this->false_multi_info[ctg][i][j]) {
							e += mi.MLClosure(local_model);
						}
						// This is <= so that, if the true structure is MFW but non-unique, there is a penalty.
						if (e <= min_energies[i]) {
							min_energies[i] = e;
							min_energy_choice[i] = static_cast<int>(j);
						}
						min_energies[i] = std::min(min_energies[i], e);
					}
				}
				double sum_fscores = 0;
				for (size_t i = 0; i < this->cts[ctg].size(); ++i) {
					if (min_energy_choice[i] == -1) {
						sum_fscores += 1;
					} else {
						sum_fscores += this->false_fscores[ctg][i][min_energy_choice[i]];
					}
				}
				double avg_error = 1.0 - (sum_fscores / this->cts[ctg].size());
				sum_sq_error += avg_error*avg_error;
			}
			return sqrt(sum_sq_error / this->cts.size());
		}, threads);

		auto it = librnary::parallel_min_element(begin(this->param_scores), end(this->param_scores), threads);
		return distance(begin(this->param_scores), it);
	}


public:

	/**
	 * @param _zero_model An energy model parameterization that gives multi-loops zero FE.
	 * @param _cts The list of CTs to use for training.
	 * @param _log Stream to log output to.
	 */
	IBFMultiLoopRMSEFscore(const ModelT &_zero_model, const VV<CTData> &_cts, std::ostream &_log)
		: IBFMultiLoop<ParamSetT, ModelT, ScorerT>(_zero_model, _cts, _log) {}
};
}


#endif //RNARK_IBF_MULTILOOP_RMSEFSCORE_HPP
