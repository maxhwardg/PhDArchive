//
// Created by max on 3/19/18.
//

#include "folders/stem_length_folder.hpp"

using namespace std;

namespace librnary {

void StemLengthFolder::SetLonelyPairs(bool v) {
	lonely_pairs = v;
}

bool StemLengthFolder::LonelyPairs() const {
	return lonely_pairs;
}

void StemLengthFolder::SetStacking(bool v) {
	stacking = v;
}

bool StemLengthFolder::Stacking() const {
	return stacking;
}

void StemLengthFolder::SetMaxTwoLoop(unsigned max_up) {
	max_twoloop_unpaired = max_up;
}

int StemLengthFolder::MaxTwoLoop() const {
	return max_twoloop_unpaired;
}

librnary::energy_t StemLengthFolder::SSScore(int i, int j) const {
	return em.Branch(i, j) + S[i][j];
}

librnary::energy_t StemLengthFolder::MLSSScore(int i, int j) const {
	return SSScore(i, j) + em.MLBranchCost();
}


void StemLengthFolder::SetModel(const StemLengthModel &_em) {
	em = _em;
}


void StemLengthFolder::TraceE(stack<TState> &s) {
	// Manage current state.
	TState curr = s.top();
	int i = curr.i;
	s.pop();

	// Base case.
	if (i == 0)
		return;

	// Maintain best decomposition.
	vector<TState> best_decomp = {TState(i - 1)}; // Unpaired 3'.
	energy_t be = E[i - 1];
	for (int k = -1; k < i; ++k) {
		int decomp = k == -1 ? 0 : E[k];
		vector<TState> decomp_state;
		if (k != -1)
			decomp_state.emplace_back(k);
		energy_t e = decomp + SSScore(k + 1, i);
		if (e < be) {
			be = e;
			best_decomp = {TState(ST, k + 1, i)};
			best_decomp.insert(best_decomp.end(), decomp_state.begin(), decomp_state.end());
		}
		if (stacking) {
			if (k + 2 < i) {
				e = decomp + SSScore(k + 2, i) + em.FiveDangle(k + 2, i);
				if (e < be) {
					be = e;
					best_decomp = {TState(ST, k + 2, i)};
					best_decomp.insert(best_decomp.end(), decomp_state.begin(), decomp_state.end());
				}
			}
			if (k + 1 < i - 1) {
				e = decomp + SSScore(k + 1, i - 1) + em.ThreeDangle(k + 1, i - 1);
				if (e < be) {
					be = e;
					best_decomp = {TState(ST, k + 1, i - 1)};
					best_decomp.insert(best_decomp.end(), decomp_state.begin(), decomp_state.end());
				}
			}
			if (k + 2 < i - 1) {
				e = decomp + SSScore(k + 2, i - 1) + em.Mismatch(k + 2, i - 1);
				if (e < be) {
					be = e;
					best_decomp = {TState(ST, k + 2, i - 1)};
					best_decomp.insert(best_decomp.end(), decomp_state.begin(), decomp_state.end());
				}
			}

			// Coaxial stack decompositions.
			for (int j = k + 1; j + 1 < i; ++j) {
				if (k + 1 < j) {
					e = decomp + em.FlushCoax(k + 1, j, j + 1, i)
						+ SSScore(k + 1, j) + SSScore(j + 1, i);
					if (e < be) {
						be = e;
						best_decomp = {TState(ST, k + 1, j), TState(ST, j + 1, i)};
						best_decomp.insert(best_decomp.end(), decomp_state.begin(), decomp_state.end());
					}
				}
				if (k + 2 < j - 1) {
					e = decomp + em.MismatchCoax(k + 2, j - 1, j + 1, i)
						+ SSScore(k + 2, j - 1) + SSScore(j + 1, i);
					if (e < be) {
						be = e;
						best_decomp = {TState(ST, k + 2, j - 1), TState(ST, j + 1, i)};
						best_decomp.insert(best_decomp.end(), decomp_state.begin(), decomp_state.end());
					}
				}
				if (j + 2 < i - 1 && k + 1 < j) {
					e = decomp + em.MismatchCoax(j + 2, i - 1, k + 1, j)
						+ SSScore(k + 1, j) + SSScore(j + 2, i - 1);
					if (e < be) {
						be = e;
						best_decomp = {TState(ST, k + 1, j), TState(ST, j + 2, i - 1)};
						best_decomp.insert(best_decomp.end(), decomp_state.begin(), decomp_state.end());
					}
				}
			}
		}

	}

	// Push best decomposition onto the stack.
	for (const auto &decomp : best_decomp) {
		s.push(decomp);
	}
}

void StemLengthFolder::TraceS(std::stack<TState> &s, Matching &m) {
	// Manage current state.
	TState curr = s.top();
	int i = curr.i, j = curr.j;
	s.pop();

	assert(i < j);

	// Maintain best decomposition.
	vector<TState> best_decomp = {TState(i - 1)};
	energy_t be = em.MaxMFE();


	// Try all helices.
	int ip = i, jp = j, sz = 1;
	energy_t fe = 0;
	while (ip < jp && ValidPair(rna[ip], rna[jp])) {
		if (fe + L[ip][jp] + em.StemLengthCost(sz) < be) {
			be = fe + L[ip][jp] + em.StemLengthCost(sz);
			best_decomp = {TState(LT, ip, jp)};
		}
		if (ip + 1 < jp - 1)
			fe += em.TwoLoop(ip, ip + 1, jp - 1, jp);
		++sz;
		++ip;
		--jp;
	}


	// Update the matching with the best stacking region.
	ip = i;
	jp = j;
	while (ip <= best_decomp[0].i) {
		m[ip] = jp;
		m[jp] = ip;
		++ip;
		--jp;
	}

	// Push best decomposition onto the stack.
	for (const auto &decomp : best_decomp) {
		s.push(decomp);
	}
}


void StemLengthFolder::TraceL(std::stack<TState> &s) {
	// Manage current state.
	TState curr = s.top();
	int i = curr.i, j = curr.j;
	s.pop();

	assert(ValidPair(rna[i], rna[j]) && (lonely_pairs || !MustBeLonelyPair(rna, i, j, em.MIN_HAIRPIN_UNPAIRED)));
	// Maintain best decomposition.
	vector<TState> best_decomp = {}; // Hairpin.
	energy_t be = em.OneLoop(i, j); // Hairpin.

	energy_t e;


	// The L (Loop) table.
	// Multi-loops.
	int init = em.MLInitCost() + em.Branch(i, j) + em.MLBranchCost();
	// No stacking interactions
	e = init + ML[2][i + 1][j - 1];
	if (e < be) {
		be = e;
		best_decomp = {TState(2, i + 1, j - 1)};
	}
	// Try stacking interactions with closing branch.
	if (stacking) {
		if (i + 2 < j - 1) { // Left dangle.
			e = init + ML[2][i + 2][j - 1] + em.ClosingThreeDangle(i, j) + em.MLUnpairedCost();
			if (e < be) {
				be = e;
				best_decomp = {TState(2, i + 2, j - 1)};
			}
		}
		if (i + 1 < j - 2) { // Right dangle.
			e = init + ML[2][i + 1][j - 2] + em.ClosingFiveDangle(i, j) + em.MLUnpairedCost();
			if (e < be) {
				be = e;
				best_decomp = {TState(2, i + 1, j - 2)};
			}
		}
		if (i + 2 < j - 2) { // Mismatch.
			e = init + ML[2][i + 2][j - 2] + em.ClosingMismatch(i, j) + em.MLUnpairedCost() * 2;
			if (e < be) {
				be = e;
				best_decomp = {TState(2, i + 2, j - 2)};
			}
		}
		// Coaxial stack.
		for (int k = i + 1; k < j; ++k) {
			// Five prime.
			// ((_)_)
			//    ^ <- k
			if (k + 1 < j - 1 && i + 1 < k) {
				e = ML[1][k + 1][j - 1] + init + em.FlushCoax(i, j, i + 1, k) + MLSSScore(i + 1, k);
				if (e < be) {
					be = e;
					best_decomp = {TState(1, k + 1, j - 1), TState(ST, i + 1, k)};
				}
			}
			// (.(_)_.)
			if (i + 2 < k && k + 1 < j - 2) {
				e = ML[1][k + 1][j - 2] + init + em.MismatchCoax(i, j, i + 2, k) + MLSSScore(i + 2, k)
					+ em.MLUnpairedCost() * 2;
				if (e < be) {
					be = e;
					best_decomp = {TState(1, k + 1, j - 2), TState(ST, i + 2, k)};
				}
			}
			// (.(_)._)
			if (i + 2 < k && k + 2 < j - 1) {
				e = ML[1][k + 2][j - 1] + init + em.MismatchCoax(i + 2, k, i, j) + MLSSScore(i + 2, k)
					+ em.MLUnpairedCost() * 2;
				if (e < be) {
					be = e;
					best_decomp = {TState(1, k + 2, j - 1), TState(ST, i + 2, k)};
				}
			}
			// Three prime.
			// (_(_))
			//   ^ <- k
			if (i + 1 < k - 1 && k < j - 1) {
				e = ML[1][i + 1][k - 1] + init + em.FlushCoax(i, j, k, j - 1) + MLSSScore(k, j - 1);
				if (e < be) {
					be = e;
					best_decomp = {TState(1, i + 1, k - 1), TState(ST, k, j - 1)};
				}
			}
			// (._(_).)
			if (k < j - 2 && i + 2 < k - 1) {
				e = ML[1][i + 2][k - 1] + init + em.MismatchCoax(i, j, k, j - 2) + MLSSScore(k, j - 2)
					+ em.MLUnpairedCost() * 2;
				if (e < be) {
					be = e;
					best_decomp = {TState(1, i + 2, k - 1), TState(ST, k, j - 2)};
				}
			}
			// (_.(_).)
			if (k < j - 2 && i + 1 < k - 2) {
				e = ML[1][i + 1][k - 2] + init + em.MismatchCoax(k, j - 2, i, j) + MLSSScore(k, j - 2)
					+ em.MLUnpairedCost() * 2;
				if (e < be) {
					be = e;
					best_decomp = {TState(1, i + 1, k - 2), TState(ST, k, j - 2)};
				}
			}

		}
	}
	// Two loops for the two-loops (bulge or internal loop).
	// Avoids stacks.
	for (int k = i + 1; k + 1 < j && (k - i - 1) <= max_twoloop_unpaired; ++k) {
		for (int l = j - 1; l > k && (j - l - 1) + (k - i - 1) <= max_twoloop_unpaired; --l) {
			// Avoid stacks.
			if ((k - i - 1) + (j - l - 1) <= 0) {
				continue;
			}
			e = S[k][l] + em.TwoLoop(i, k, l, j);
			if (e < be) {
				be = e;
				best_decomp = {TState(ST, k, l)};
			}
		}
	}

	// Push best decomposition onto the stack.
	for (const auto &decomp : best_decomp) {
		s.push(decomp);
	}
}

void StemLengthFolder::TraceCx(std::stack<TState> &s) {
	// Manage current state.
	TState curr = s.top();
	int i = curr.i, j = curr.j;
	s.pop();

	assert(i < j);
	assert(stacking);

	// Maintain best decomposition.
	vector<TState> best_decomp = {};
	energy_t be = em.MaxMFE(), e;

	for (int k = i + 1; k + 1 < j; ++k) {
		e = em.FlushCoax(i, k, k + 1, j) + MLSSScore(i, k) + MLSSScore(k + 1, j);
		if (e < be) {
			be = e;
			best_decomp = {TState(ST, i, k), TState(ST, k + 1, j)};
		}
		if (i + 1 < k - 1) {
			e = em.MismatchCoax(i + 1, k - 1, k + 1, j) + MLSSScore(i + 1, k - 1) + MLSSScore(k + 1, j)
				+ em.MLUnpairedCost() * 2;
			if (e < be) {
				be = e;
				best_decomp = {TState(ST, i + 1, k - 1), TState(ST, k + 1, j)};
			}
		}
		if (k + 2 < j - 1) {
			e = em.MismatchCoax(k + 2, j - 1, i, k) + MLSSScore(i, k) + MLSSScore(k + 2, j - 1)
				+ em.MLUnpairedCost() * 2;
			if (e < be) {
				be = e;
				best_decomp = {TState(ST, i, k), TState(ST, k + 2, j - 1)};
			}
		}
	}

	// Push best decomposition onto the stack.
	for (const auto &decomp : best_decomp) {
		s.push(decomp);
	}
}

void StemLengthFolder::TraceML(std::stack<TState> &s) {
	// Manage current state.
	TState curr = s.top();
	int i = curr.i, j = curr.j, b = curr.b;
	s.pop();

	if (i == j) {
		assert(b == 0);
		return;
	}

	assert(i < j);
	assert(b >= 0 && b <= 2);

	// Maintain best decomposition.
	vector<TState> best_decomp = {TState(b, i, j - 1)}; // Unpaired on the right.
	energy_t be = ML[b][i][j - 1] + em.MLUnpairedCost(), e;

	if (b < 2) { // End on branch cases.
		e = MLSSScore(i, j);
		if (e < be) {
			be = e;
			best_decomp = {TState(ST, i, j)};
		}
		if (stacking) {
			if (i + 1 < j) {
				e = MLSSScore(i + 1, j) + em.FiveDangle(i + 1, j) + em.MLUnpairedCost();
				if (e < be) {
					be = e;
					best_decomp = {TState(ST, i + 1, j)};
				}
				e = MLSSScore(i, j - 1) + em.ThreeDangle(i, j - 1) + em.MLUnpairedCost();
				if (e < be) {
					be = e;
					best_decomp = {TState(ST, i, j - 1)};
				}
			}
			if (i + 1 < j - 1) {
				e = MLSSScore(i + 1, j - 1) + em.Mismatch(i + 1, j - 1) + em.MLUnpairedCost() * 2;
				if (e < be) {
					be = e;
					best_decomp = {TState(ST, i + 1, j - 1)};
				}
			}
		}
	}
	// End on coaxial stack.
	if (stacking) {
		e = Cx[i][j];
		if (e < be) {
			be = e;
			best_decomp = {TState(CxT, i, j)};
		}
	}
	// bprime is the number of branches required after one has been placed.
	int bprime = max(0, b - 1);
	for (int k = i; k + 2 <= j; ++k) { // Try all decompositions into 5' ML fragment and 3' branch.
		e = ML[bprime][i][k] + MLSSScore(k + 1, j);
		if (e < be) {
			be = e;
			best_decomp = {TState(bprime, i, k), TState(ST, k + 1, j)};
		}
		// From here on is stacking.
		if (stacking) {
			if (k + 2 < j) {
				e = ML[bprime][i][k] + MLSSScore(k + 2, j) + em.FiveDangle(k + 2, j) + em.MLUnpairedCost();
				if (e < be) {
					be = e;
					best_decomp = {TState(bprime, i, k), TState(ST, k + 2, j)};
				}
			}
			if (k + 1 < j - 1) {
				e = ML[bprime][i][k] + MLSSScore(k + 1, j - 1) + em.ThreeDangle(k + 1, j - 1) + em.MLUnpairedCost();
				if (e < be) {
					be = e;
					best_decomp = {TState(bprime, i, k), TState(ST, k + 1, j - 1)};
				}
			}
			if (k + 2 < j - 1) {
				e = ML[bprime][i][k] + MLSSScore(k + 2, j - 1) + em.Mismatch(k + 2, j - 1) + em.MLUnpairedCost() * 2;
				if (e < be) {
					be = e;
					best_decomp = {TState(bprime, i, k), TState(ST, k + 2, j - 1)};
				}
			}
			// Coaxial stack decomposition.
			e = ML[0][i][k] + Cx[k + 1][j];
			if (e < be) {
				be = e;
				best_decomp = {TState(0, i, k), TState(CxT, k + 1, j)};
			}
		}
	}

	// Push best decomposition onto the stack.
	for (const auto &decomp : best_decomp) {
		s.push(decomp);
	}
}


librnary::Matching StemLengthFolder::Traceback() {
	const int N = static_cast<int>(rna.size());
	Matching m = EmptyMatching(static_cast<unsigned>(rna.size()));
	if (N == 0) {
		return m;
	}
	stack<TState> s;
	s.push(TState(N - 1));
	while (!s.empty()) {
		if (s.top().t == ET)
			TraceE(s);
		else if (s.top().t == ST) {
			TraceS(s, m);
		} else if (s.top().t == LT) {
			TraceL(s);
		} else if (s.top().t == MLT) {
			TraceML(s);
		} else // s.top().t == Cx
			TraceCx(s);
	}
	return m;
}

librnary::energy_t StemLengthFolder::Fold(const PrimeStructure &_rna) {
	// Load the RNA into the energy model.
	em.SetRNA(_rna);
	// Save the RNA.
	this->rna = _rna;

	// The length of the RNA.
	const unsigned long RSZ = rna.size();
	// N is a convenience variable.
	// I have it as well as RSZ to avoid type warnings.
	// At least forcing this cast will catch any overflow errors, potentially.
	const int N = static_cast<int>(RSZ);

	if (N == 0) {
		return 0;
	}

	// Clear the DP tables.
	S.clear();
	L.clear();
	ML.clear();
	Cx.clear();
	E.clear();

	// Initialize the DP tables.
	S.resize(RSZ, VE(RSZ, em.MaxMFE()));
	Cx = L = S;
	ML.resize(3, VVE(RSZ, VE(RSZ, em.MaxMFE())));
	E.resize(RSZ, 0);

	// Special base for for ML. End on a single unpaired.
	for (int i = 0; i < N; ++i)
		ML[0][i][i] = em.MLUnpairedCost();

	for (int i = N - 2; i >= 0; --i) { // i is 5' nucleotide.
		for (int j = i + 1; j < N; ++j) { // j is 3' nucleotide.
			if (ValidPair(rna[i], rna[j]) &&
				(lonely_pairs || !MustBeLonelyPair(rna, i, j, em.MIN_HAIRPIN_UNPAIRED))) {
				// The L (Loop) table.
				energy_t best = em.OneLoop(i, j); // Hairpin.
				// Multi-loops.
				int init = em.MLInitCost() + em.Branch(i, j) + em.MLBranchCost();
				// No stacking interactions
				best = min(best, init + ML[2][i + 1][j - 1]);
				// Try stacking interactions with closing branch.
				if (stacking) {
					if (i + 2 < j - 1) // Left dangle.
						best = min(best,
								   init + ML[2][i + 2][j - 1] + em.ClosingThreeDangle(i, j) + em.MLUnpairedCost());
					if (i + 1 < j - 2) // Right dangle.
						best = min(best,
								   init + ML[2][i + 1][j - 2] + em.ClosingFiveDangle(i, j) + em.MLUnpairedCost());
					if (i + 2 < j - 2) // Mismatch.
						best = min(best,
								   init + ML[2][i + 2][j - 2] + em.ClosingMismatch(i, j) + em.MLUnpairedCost() * 2);
					// Coaxial stack.
					for (int k = i + 1; k < j; ++k) {
						// Five prime.
						// ((_)_)
						//    ^ <- k
						if (k + 1 < j - 1 && i + 1 < k)
							best = min(best,
									   ML[1][k + 1][j - 1] + init + em.FlushCoax(i, j, i + 1, k)
										   + MLSSScore(i + 1, k));
						// (.(_)_.)
						if (i + 2 < k && k + 1 < j - 2)
							best = min(best,
									   ML[1][k + 1][j - 2] + init + em.MismatchCoax(i, j, i + 2, k)
										   + MLSSScore(i + 2, k) + em.MLUnpairedCost() * 2);
						// (.(_)._)
						if (i + 2 < k && k + 2 < j - 1)
							best = min(best,
									   ML[1][k + 2][j - 1] + init + em.MismatchCoax(i + 2, k, i, j)
										   + MLSSScore(i + 2, k) + em.MLUnpairedCost() * 2);
						// Three prime.
						// (_(_))
						//   ^ <- k
						if (i + 1 < k - 1 && k < j - 1)
							best = min(best,
									   ML[1][i + 1][k - 1] + init + em.FlushCoax(i, j, k, j - 1)
										   + MLSSScore(k, j - 1));
						// (._(_).)
						if (k < j - 2 && i + 2 < k - 1)
							best = min(best,
									   ML[1][i + 2][k - 1] + init + em.MismatchCoax(i, j, k, j - 2)
										   + MLSSScore(k, j - 2) + em.MLUnpairedCost() * 2);
						// (_.(_).)
						if (k < j - 2 && i + 1 < k - 2)
							best = min(best,
									   ML[1][i + 1][k - 2] + init + em.MismatchCoax(k, j - 2, i, j)
										   + MLSSScore(k, j - 2) + em.MLUnpairedCost() * 2);
					}
				}
				// Two loops for the two-loops (bulge or internal loop).
				// Avoids stacks,
				for (int k = i + 1; k + 1 < j && (k - i - 1) <= max_twoloop_unpaired; ++k) {
					for (int l = j - 1; l > k && (j - l - 1) + (k - i - 1) <= max_twoloop_unpaired; --l) {
						// Avoid stacks.
						if ((k - i - 1) + (j - l - 1) <= 0) {
							continue;
						}
						best = min(best, S[k][l] + em.TwoLoop(i, k, l, j));
					}
				}
				L[i][j] = best;

				// Now the S (Stacking) table.
				best = em.MaxMFE();
				// Try all helices.
				int ip = i, jp = j, sz = 1;
				energy_t fe = 0;
				while (ip < jp && ValidPair(rna[ip], rna[jp])) {
					best = min(best, fe + L[ip][jp] + em.StemLengthCost(sz));
					if (ip + 1 < jp - 1)
						fe += em.TwoLoop(ip, ip + 1, jp - 1, jp);
					++sz;
					++ip;
					--jp;
				}
				S[i][j] = best;
			}

			// Fill the coaxial stack table if stacking is enabled.
			energy_t best = em.MaxMFE();
			for (int k = i + 1; k + 1 < j && stacking; ++k) {
				best = min(best, em.FlushCoax(i, k, k + 1, j) + MLSSScore(i, k) + MLSSScore(k + 1, j));
				if (i + 1 < k - 1)
					best = min(best,
							   em.MismatchCoax(i + 1, k - 1, k + 1, j) + MLSSScore(i + 1, k - 1) + MLSSScore(k + 1, j)
								   + em.MLUnpairedCost() * 2);
				if (k + 2 < j - 1)
					best = min(best,
							   em.MismatchCoax(k + 2, j - 1, i, k) + MLSSScore(i, k) + MLSSScore(k + 2, j - 1)
								   + em.MLUnpairedCost() * 2);
			}
			Cx[i][j] = best;

			for (int b = 0; b < 3; ++b) { // b is the branches needed for valid ML.
				best = ML[b][i][j - 1] + em.MLUnpairedCost();
				if (b < 2) { // End on branch cases.
					best = min(best, MLSSScore(i, j));
					if (stacking) {
						if (i + 1 < j) {
							best = min(best, MLSSScore(i + 1, j) + em.FiveDangle(i + 1, j) + em.MLUnpairedCost());
							best = min(best, MLSSScore(i, j - 1) + em.ThreeDangle(i, j - 1) + em.MLUnpairedCost());
						}
						if (i + 1 < j - 1) {
							best = min(best,
									   MLSSScore(i + 1, j - 1) + em.Mismatch(i + 1, j - 1) + em.MLUnpairedCost() * 2);
						}
					}
				}
				// End on coaxial stack.
				if (stacking) {
					best = min(best, Cx[i][j]);
				}
				// bprime is the number of branches required after one has been placed.
				int bprime = max(0, b - 1);
				for (int k = i; k + 2 <= j; ++k) { // Try all decompositions into 5' ML fragment and 3' branch.
					best = min(best, ML[bprime][i][k] + MLSSScore(k + 1, j));
					// From here on is stacking.
					if (stacking) {
						if (k + 2 < j)
							best = min(best,
									   ML[bprime][i][k] + MLSSScore(k + 2, j) + em.FiveDangle(k + 2, j)
										   + em.MLUnpairedCost());
						if (k + 1 < j - 1)
							best = min(best,
									   ML[bprime][i][k] + MLSSScore(k + 1, j - 1) + em.ThreeDangle(k + 1, j - 1)
										   + em.MLUnpairedCost());
						if (k + 2 < j - 1)
							best = min(best,
									   ML[bprime][i][k] + MLSSScore(k + 2, j - 1) + em.Mismatch(k + 2, j - 1)
										   + em.MLUnpairedCost() * 2);
						// Coaxial stack decomposition.
						best = min(best, ML[0][i][k] + Cx[k + 1][j]);
					}
				}
				ML[b][i][j] = best;
			}
		}
	}


	for (int i = 1; i < N; ++i) {
		energy_t best = E[i - 1];
		for (int k = -1; k < i; ++k) {
			energy_t decomp = k == -1 ? 0 : E[k];
			best = min(best, decomp + SSScore(k + 1, i));
			if (stacking) {
				if (k + 2 < i)
					best = min(best, decomp + SSScore(k + 2, i) + em.FiveDangle(k + 2, i));
				if (k + 1 < i - 1)
					best = min(best, decomp + SSScore(k + 1, i - 1) + em.ThreeDangle(k + 1, i - 1));
				if (k + 2 < i - 1)
					best = min(best, decomp + SSScore(k + 2, i - 1) + em.Mismatch(k + 2, i - 1));
				// Coaxial stack decompositions.
				for (int j = k + 1; j + 1 < i; ++j) {
					if (k + 1 < j)
						best = min(best, decomp + em.FlushCoax(k + 1, j, j + 1, i)
							+ SSScore(k + 1, j) + SSScore(j + 1, i));
					if (k + 2 < j - 1)
						best = min(best, decomp + em.MismatchCoax(k + 2, j - 1, j + 1, i)
							+ SSScore(k + 2, j - 1) + SSScore(j + 1, i));
					if (k + 1 < j && j + 2 < i - 1)
						best = min(best, decomp + em.MismatchCoax(j + 2, i - 1, k + 1, j)
							+ SSScore(k + 1, j) + SSScore(j + 2, i - 1));
				}
			}
		}
		E[i] = best;
	}

	return E[N - 1];
}

}