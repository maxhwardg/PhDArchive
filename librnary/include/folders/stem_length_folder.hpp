//
// Created by max on 3/19/18.
//

#ifndef RNARK_STEM_LENGTH_FOLDER_HPP
#define RNARK_STEM_LENGTH_FOLDER_HPP

#include <stack>

#include "vector_types.hpp"
#include "primary_structure.hpp"
#include "models/stem_length_model.hpp"


namespace librnary {
/**
 * This class does MFE RNA folding using the typical Turner-Mathews nearest neighbour rules. However, these have been
 * modified slightly so that stems have a free energy change cost based on their length. See StemLengthModel for more
 * details.
 */
class StemLengthFolder {
protected:
	PrimeStructure rna;
	StemLengthModel em;

	/**
	 * The 'Stem' table. S[i][j] is the optimal substructure closed by a pair i,j such that i,j are the start of a stem.
	 */
	VVE S;

	/**
	 * The 'Loop' table. L[i][j] is the optimal substructure closed by the pair i,j such that i,j close a loop.
	 * Note that a stacking region is not a loop.
	 */
	VVE L;

	/**
	 * The 'Multi-Loop' table. ML[b][i][j] is the optimal part of the multi-loop that definitely has
	 * at least b branches.
	 */
	VVVE ML;

	/**
	 * The 'Coaxial stack' table. Cx[i][j] is the optimal multi-loop coaxial stack such that one branch starts at i,
	 * and the other branch ends at j. Note that it is assumed that these branches are in a multi-loop. This table will
	 * count the unpaired cost.
	 */
	VVE Cx;

	/**
	 * The 'External loop' table. E[i] is the optimal external loop fragment 0..i.
	 */
	VE E;

	/**
	 * Represents a particular table during a traceback.
	 */
	enum Table {
		ET, ST, LT, MLT, CxT
	};

	/**
	 * Trace state during traceback. Represents a location in a table.
	 */
	struct TState {
		Table t;
		int b, i, j;

		TState(int _i)
			: t(ET), i(_i) {};

		TState(Table _t, int _i, int _j)
			: t(_t), i(_i), j(_j) {};

		TState(int _b, int _i, int _j)
			: t(MLT), b(_b), i(_i), j(_j) {};
	};

	/**
	 * Traceback through a place in the E table. Assumes this place/state is at the top of s.
	 * Updates s accordingly with the decomposed states.
	 * @param s Current open states in the traceback in stack form.
	 */
	void TraceE(std::stack<TState> &s);

	/**
	 * Traceback through a place in the S table. Assumes this place/state is at the top of s.
	 * Updates s accordingly with the decomposed states. Also updates the matching with pairs in the stem.
	 * @param s Current open states in the traceback in stack form.
	 * @param m The current partially traced structure. Will modify this by adding pairs to it.
	 */
	void TraceS(std::stack<TState> &s, Matching &m);


	/**
	 * Traceback through a place in the L table. Assumes this place/state is at the top of s.
	 * Updates s accordingly with the decomposed states.
	 * @param s Current open states in the traceback in stack form.
	 */
	void TraceL(std::stack<TState> &s);


	/**
	 * Traceback through a place in the Cx table. Assumes this place/state is at the top of s.
	 * Updates s accordingly with the decomposed states.
	 * @param s Current open states in the traceback in stack form.
	 */
	void TraceCx(std::stack<TState> &s);

	/**
	 * Traceback through a place in the ML table. Assumes this place/state is at the top of s.
	 * Updates s accordingly with the decomposed states.
	 * @param s Current open states in the traceback in stack form.
	 */
	void TraceML(std::stack<TState> &s);

	/**
	 * This flag toggles whether lonely pairs are allowed.
	 * Lonely pairs is a bit weird, since a lonely pair in the sense of this flag isn't a length-1 stem. It is actually
	 * a pair that MUST be a length-1 stem due to the adjacent nucleotides.
	 */
	bool lonely_pairs = true;

	/// This flag toggles whether stacking interactions (dangles, terminal mismatch, and coaxial stacking) are used.
	bool stacking = true;

	/// The maximum number of unpaired in a bulge or internal loop.
	int max_twoloop_unpaired = std::numeric_limits<int>::max() / 3;

	/**
	 * The optimal sub-surface score of the structure closed by (i,j). Designed for external-loop branches.
	 * Accounts for AU/GU penalty.
	 */
	energy_t SSScore(int i, int j) const;

	/**
	 * The optimal sub-surface score of the structure closed by (i,j) in a multi-loop.
	 * Accounts for AU/GU penalty, and multi-loop branch cost.
	 */
	energy_t MLSSScore(int i, int j) const;

public:

	/// Set the max unpaired nucleotides in a bulge/internal loop.
	void SetMaxTwoLoop(unsigned max_up);

	/// Get the max unpaired nucleotides in a bulge/internal loop.
	int MaxTwoLoop() const;

	void SetStacking(bool v);

	bool Stacking() const;

	void SetLonelyPairs(bool v);

	bool LonelyPairs() const;

	void SetModel(const StemLengthModel &_em);

	energy_t Fold(const PrimeStructure &_rna);

	/**
	 * Trace back through the DP tables and produce a MFE secondary structure.
	 * @return The traced secondary structure.
	 */
	Matching Traceback();

	StemLengthFolder(const StemLengthModel &_em)
		: em(_em) {}

};

}

#endif //RNARK_STEM_LENGTH_FOLDER_HPP
