//
// Created by max on 3/19/18.
//

#include "models/stem_length_model.hpp"

namespace librnary {

energy_t StemLengthModel::StemLengthCost(int length) const {
	assert(length > 0);
	if (static_cast<int>(stem_length_costs.size()) < length)
		return stem_length_costs.back();
	return stem_length_costs[length - 1];
}

void StemLengthModel::SetLengthCosts(const std::vector<energy_t> &length_costs) {
	stem_length_costs = length_costs;
}

std::vector<energy_t> StemLengthModel::StemLengthCosts() const {
	return stem_length_costs;
}

}