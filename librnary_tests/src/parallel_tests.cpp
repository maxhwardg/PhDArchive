//
// Created by max on 8/31/17.
//

#include <gtest/gtest.h>

#include "parallel.hpp"
#include "random.hpp"

using namespace std;

TEST(Parallel, SmallVectorsTransform) {
	for (size_t TC = 0; TC < 100; ++TC) {
		vector<long> A(17);
		auto re = librnary::RandomEngineForTests();
		for (size_t i = 0; i < A.size(); ++i) {
			A[i] = re();
		}
		long num = re();
		vector<long> res(A.size()), res2(A.size());
		std::transform(A.begin(), A.end(), res.begin(), [&](long e) {
			return e+num;
		});
		librnary::parallel_transform(A, res2, [&](long e) {
			return e+num;
		});
		EXPECT_EQ(res, res2);
	}
}

TEST(Parallel, LargeVectorsTransform) {
	for (size_t TC = 0; TC < 100; ++TC) {
		vector<long> A(523);
		auto re = librnary::RandomEngineForTests();
		for (size_t i = 0; i < A.size(); ++i) {
			A[i] = re();
		}
		long num = re();
		vector<long> res(A.size()), res2(A.size());
		std::transform(A.begin(), A.end(), res.begin(), [&](long e) {
			return e+num;
		});
		librnary::parallel_transform(A, res2, [&](long e) {
			return e+num;
		});
		EXPECT_EQ(res, res2);
	}
}

TEST(Parallel, SmallVectorsMinMax) {
	for (size_t TC = 0; TC < 100; ++TC) {
		vector<long> A(17);
		auto re = librnary::RandomEngineForTests();
		for (size_t i = 0; i < A.size(); ++i) {
			A[i] = re();
		}
		EXPECT_EQ(std::min_element(A.begin(), A.end()), librnary::parallel_min_element(A.begin(), A.end()));
		EXPECT_EQ(std::max_element(A.begin(), A.end()), librnary::parallel_max_element(A.begin(), A.end()));
	}
}

TEST(Parallel, LargeVectorsMinMax) {
	vector<long> A(523);
	for (size_t TC = 0; TC < 100; ++TC) {
		auto re = librnary::RandomEngineForTests();
		for (size_t i = 0; i < A.size(); ++i) {
			A[i] = re();
		}
		EXPECT_EQ(std::max_element(A.begin(), A.end()), librnary::parallel_max_element(A.begin(), A.end()));
		EXPECT_EQ(std::min_element(A.begin(), A.end()), librnary::parallel_min_element(A.begin(), A.end()));
	}
}