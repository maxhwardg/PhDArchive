//
// Created by max on 10/23/17.
//

#ifndef RNARK_IBG_MULTILOOP_ANDRONESCU_HPP
#define RNARK_IBG_MULTILOOP_ANDRONESCU_HPP

#include "IBF_multiloop.hpp"
namespace librnary {
/**
 * A extenstion of the typical IBF trainer that emulates the Andronescu algorithm.
 * The the Andronescu algorithm appears to minimize the sum of squared errors for each RNA.
 * @tparam ParamSetT Parameter set type.
 * @tparam ModelT Energy model type.
 * @tparam ScorerT Energy model socrer type.
 */
template<typename ParamSetT, typename ModelT, typename ScorerT, typename FolderT>
class IBFMultiLoopAndronescu: public IBFMultiLoop<ParamSetT, ModelT, ScorerT, FolderT> {

	long FindBestParams(const std::vector<ParamSetT> &params) override {
		librnary::parallel_transform(params, this->param_scores, [=](const ParamSetT &pset) {
			auto local_model = this->zero_model;
			pset.LoadInto(local_model);

			double sum_sq_errors = 0;

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
						// This is <= so that, if the true structure is MFE but non-unique, there is a penalty.
						if (e <= min_energies[i]) {
							min_energies[i] = e;
							min_energy_choice[i] = static_cast<int>(j);
						}
					}
				}
				for (size_t i = 0; i < this->cts[ctg].size(); ++i) {
					if (min_energy_choice[i] != -1) {
						int error = abs(true_energies[i] - min_energies[i]);
						sum_sq_errors += error * error;
					}
				}
			}
			return sum_sq_errors;
		}, this->threads);

		auto it = librnary::parallel_min_element(begin(this->param_scores), end(this->param_scores), this->threads);
		return distance(begin(this->param_scores), it);
	}
public:

	/**
	 * @param _zero_model An energy model parameterization that gives multi-loops zero FE.
	 * @param _cts The list of CTs to use for training.
	 * @param _log Stream to log output to.
	 */
	IBFMultiLoopAndronescu(const ModelT &_zero_model, VV<CTData> _cts, std::ostream &_log)
		: IBFMultiLoop<ParamSetT, ModelT, ScorerT, FolderT>(_zero_model, _cts, _log) {}
};

}

#endif //RNARK_IBG_MULTILOOP_ANDRONESCU_HPP
