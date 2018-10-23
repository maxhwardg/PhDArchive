//
// Created by max on 3/19/18.
//

#ifndef RNARK_STEM_LENGTH_MODEL_HPP
#define RNARK_STEM_LENGTH_MODEL_HPP

#include "nn_affine_model.hpp"

namespace librnary {

/**
 * A model that extends the NNAffineModel and adds a cost to stems based on their length. This can be thought of as a
 * generalization of the no lonely pairs heuristic. A cost can be given to stems with length 1 (a lonely pair), and a
 * different cost can be given to stems of other sizes too.
 */
class StemLengthModel: public NNAffineModel {
protected:
	std::vector<energy_t> stem_length_costs = {0};
public:
	/**
	 * Computes the free energy change cost of a stem.
	 * @param length Length of the stem in pairs.
	 * @return Free energy change.
	 */
	energy_t StemLengthCost(int length) const;
	/**
	 * Assigns a new stem length costs vector. Index 0 is the cost for a stem of length 1, index 1 is the cost for a
	 * stem of length 2, and so on. Any stem whose length would be an invalid index will use the last valid index's
	 * cost.
	 * @param length_costs The new stem length costs vector.
	 */
	void SetLengthCosts(const std::vector<energy_t> &length_costs);
	/**
	 * Gets the stem length costs vector. See SetLengthCosts.
	 * @return The stem length costs vector.
	 */
	std::vector<energy_t> StemLengthCosts() const;
	StemLengthModel(const std::string &data_path, const PrimeStructure &_rna)
		: NNAffineModel(data_path, _rna) {}
	StemLengthModel(const std::string &data_path)
		: NNAffineModel(data_path) {}
};

}

#endif //RNARK_STEM_LENGTH_MODEL_HPP
