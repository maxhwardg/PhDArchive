//
// Created by max on 5/30/16.
//

#include <chrono>
#include <iostream>

#include "random.hpp"

using namespace std;

std::default_random_engine librnary::RandomEngineForTests() {
	// The seed is hours since epoch. This provides a nice slowly changing seed.
	long seed = chrono::duration_cast<chrono::hours>(chrono::system_clock::now().time_since_epoch()).count();
	// If a test fails, the last seed before it should reproduce the problem.
	clog << "Random Engine Seed: " << seed << endl;
	default_random_engine re(seed);
	return re;
}

