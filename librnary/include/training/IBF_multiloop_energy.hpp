//
// Created by max on 7/27/17.
//

#ifndef RNARK_IBF_MULTILOOP_ENERGY_HPP
#define RNARK_IBF_MULTILOOP_ENERGY_HPP

#include "training/IBF_multiloop.hpp"
#include "parallel.hpp"

namespace librnary {

/**
 * A modification of IBFMultiLoop that uses RMSE of energy differences between the true structure and the MFE structure
 * as a guide to parameter fitness.
 * Deals with sets by taking the RMS of RMSEs for each set.
 */
template<typename ParamSetT, typename ModelT, typename ScorerT>
class IBFMultiLoopEnergy : public IBFMultiLoop<ParamSetT, ModelT, ScorerT> {
private:

protected:
	virtual long FindBestParams(const std::vector<ParamSetT> &params, const ModelT &model) {
		librnary::parallel_transform(params, this->param_scores, [=](const ParamSetT &pset) {
			auto local_model = model;
			pset.LoadInto(local_model);

			double sum_sq_rmses = 0;

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
				double sum_sq_error_energy = 0;
				for (size_t i = 0; i < this->cts[ctg].size(); ++i) {
					if (min_energies[i] == 0)
						min_energies[i] = 1;
					double scaled = (true_energies[i] - min_energies[i]) / static_cast<double>(true_energies[i]);
					sum_sq_error_energy += scaled * scaled;
				}
				sum_sq_rmses += sqrt(sum_sq_error_energy/this->cts[ctg].size());
			}
			return sqrt(sum_sq_rmses / this->cts.size());
		}, threads);

		auto it = librnary::parallel_min_element(begin(this->param_scores), end(this->param_scores), threads);
		return distance(begin(this->param_scores), it);
	}


public:

	IBFMultiLoopEnergy(ModelT _zero_model, VV<CTData> _cts, std::ostream &_log)
		: IBFMultiLoop<ParamSetT, ModelT, ScorerT>(_zero_model, _cts, _log) {}
};
}

#endif //RNARK_IBF_MULTILOOP_ENERGY_HPP
