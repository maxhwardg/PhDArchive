//
// Created by max on 5/30/16.
//
/*
 * Contains functions for random numbers used by librnary tests.
 */

#ifndef RNARK_RANDOM_HPP
#define RNARK_RANDOM_HPP

#include <random>

namespace librnary {
/*
 * This function returns a default random engine for use in testing.
 * It is good to use this function, as it changes the seed, thus leading to more coverage over time.
 */
std::default_random_engine RandomEngineForTests();
}

#endif //RNARK_RANDOM_HPP
