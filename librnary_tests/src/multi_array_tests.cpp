//
// Created by max on 6/23/16.
//


#include <gtest/gtest.h>
#include <multi_array.hpp>
#include <vector_types.hpp>
#include "random.hpp"

using namespace std;

TEST(MultiArray, MultiArrayBasicTest) {
	librnary::MultiArray<int, 2> arr({3, 3}, 0);
	int v = arr[{0, 0}];
	EXPECT_EQ(v, 0);
	arr[{0, 0}] = 1;
	arr[{0, 2}] = 3;
	v = arr[{0, 0}];
	EXPECT_EQ(v, 1);
	v = arr[{0, 2}];
	EXPECT_EQ(v, 3);
	arr[{1, 1}] = 5;
	arr[{1, 2}] = 6;
	v = arr[{0, 2}];
	EXPECT_EQ(v, 3);
	v = arr[{0, 0}];
	EXPECT_EQ(v, 1);
	v = arr[{1, 2}];
	EXPECT_EQ(v, 6);
}

TEST(MultiArray, Array2DBasicTest) {
	librnary::Array2D<int> arr(3, 3, 0);
	EXPECT_EQ(arr[0][0], 0);
	arr[0][0] = 1;
	arr[0][2] = 3;
	EXPECT_EQ(arr[0][0], 1);
	EXPECT_EQ(arr[0][2], 3);
	arr[1][1] = 5;
	arr[1][2] = 6;
	EXPECT_EQ(arr[0][2], 3);
	EXPECT_EQ(arr[0][0], 1);
	EXPECT_EQ(arr[1][2], 6);
	EXPECT_DEBUG_DEATH(arr[6][1], "");
	EXPECT_DEBUG_DEATH(arr[0][3], "");
}

TEST(MultiArray, MultiArrayFuzzAgainstVector3D) {
	const int A = 23, B = 33, C = 51, CASES = 50000;
	librnary::VVVI vec(A, librnary::VVI(B, librnary::VI(C, 0)));
	librnary::MultiArray<int, 3> arr({A, B, C}, 0);
	auto re = librnary::RandomEngineForTests();
	for (int tc = 0; tc < CASES; ++tc) {
		int a = re() % A, b = re() % B, c = re() % C, v = re() % 100 - 50;
		EXPECT_EQ(vec[a][b][c], arr(a, b, c));
		vec[a][b][c] = v;
		arr(a, b, c) = v;
	}
}

TEST(MultiArray, Array3DFuzzAgainstVector3D) {
	const int A = 23, B = 33, C = 51, CASES = 50000;
	librnary::VVVI vec(A, librnary::VVI(B, librnary::VI(C, 0)));
	librnary::Array3D<int> arr(A, B, C, 0);
	auto re = librnary::RandomEngineForTests();
	for (int tc = 0; tc < CASES; ++tc) {
		int a = re() % A, b = re() % B, c = re() % C, v = re() % 100 - 50;
		EXPECT_EQ(vec[a][b][c], arr[a][b][c]);
		vec[a][b][c] = v;
		arr[a][b][c] = v;
	}
}

TEST(MultiArray, MultiArrayFuzzAgainstVector4D) {
	const int A = 11, B = 12, C = 33, D = 39, CASES = 50000;
	librnary::_4DVI vec(A, librnary::VVVI(B, librnary::VVI(C, librnary::VI(D))));
	librnary::MultiArray<int, 4> arr({A, B, C, D}, 0);
	auto re = librnary::RandomEngineForTests();
	for (int tc = 0; tc < CASES; ++tc) {
		int a = re() % A, b = re() % B, c = re() % C, d = re() % D, v = re() % 100 - 50;
		EXPECT_EQ(vec[a][b][c][d], arr(a, b, c, d));
		vec[a][b][c][d] = v;
		arr(a, b, c, d) = v;
	}
//	cout << vec[0][0][0][0] << endl;
//	cout << arr(0,0,0,0) << endl;
}

TEST(MultiArray, Array4DFuzzAgainstVector4D) {
	const int A = 8, B = 12, C = 33, D = 39, CASES = 50000;
	librnary::_4DVI vec(A, librnary::VVVI(B, librnary::VVI(C, librnary::VI(D))));
	librnary::Array4D<int> arr(A, B, C, D, 0);
	auto re = librnary::RandomEngineForTests();
	for (int tc = 0; tc < CASES; ++tc) {
		int a = re() % A, b = re() % B, c = re() % C, d = re() % D, v = re() % 100 - 50;
		EXPECT_EQ(vec[a][b][c][d], arr[a][b][c][d]);
		vec[a][b][c][d] = v;
		arr[a][b][c][d] = v;
	}
//	cout << vec[0][0][0][0] << endl;
//	cout << arr[0][0][0][0] << endl;
}

TEST(MultiArray, MultiArrayIterateAgainstVector4D) {
	const int A = 15, B = 30, C = 50, D = 50;
	librnary::_4DVI vec(A, librnary::VVVI(B, librnary::VVI(C, librnary::VI(D))));
	librnary::MultiArray<int, 4> arr({A, B, C, D}, 0);
	for (int c = C - 1; c >= 0; --c) {
		for (int d = c; d < D; ++d) {
			for (int a = 0; a < A; ++a) {
				for (int b = 0; b < B; ++b) {
					for (int i = c + 1; i < d; ++i) {
						vec[a][b][i - 1][d] += 2;
						vec[a][b][c][i - 1] += 2;
						arr(a, b, i - 1, d) += 2;
						arr(a, b, c, i - 1) += 2;
					}
					EXPECT_EQ(arr(a, b, c, d), vec[a][b][c][d]);
				}
			}
		}
	}
}

TEST(MultiArray, Array4DIterateAgainstVector4D) {
	const int A = 15, B = 30, C = 50, D = 50;
	librnary::_4DVI vec(A, librnary::VVVI(B, librnary::VVI(C, librnary::VI(D))));
	librnary::Array4D<int> arr(A, B, C, D, 0);
	for (int c = C - 1; c >= 0; --c) {
		for (int d = c; d < D; ++d) {
			for (int a = 0; a < A; ++a) {
				for (int b = 0; b < B; ++b) {
					for (int i = c + 1; i < d; ++i) {
						vec[a][b][i - 1][d] += 2;
						vec[a][b][c][i - 1] += 2;
						arr[a][b][i - 1][d] += 2;
						arr[a][b][c][i - 1] += 2;
					}
					EXPECT_EQ(arr[a][b][c][d], vec[a][b][c][d]);
				}
			}
		}
	}
}

TEST(MultiArray, CopyTest1D) {
	librnary::Array1D<int> arr(10, 8);
	arr[2] = 3;
	auto test = arr;
	EXPECT_EQ(test[2], 3);
	EXPECT_EQ(test[9], 8);
}

TEST(MultiArray, CopyTest2D) {
	librnary::Array2D<int> arr(10, 11, 8);
	arr[2][8] = 3;
	auto test = arr;
	EXPECT_EQ(test[2][8], 3);
	EXPECT_EQ(test[9][7], 8);
}

TEST(MultiArray, CopyTest3D) {
	librnary::Array3D<int> arr(10, 11, 9, 8);
	arr[2][8][1] = 3;
	auto test = arr;
	EXPECT_EQ(test[2][8][1], 3);
	EXPECT_EQ(test[9][7][8], 8);
}

TEST(MultiArray, CopyTest4D) {
	librnary::Array4D<int> arr(10, 11, 9, 13, 8);
	arr[2][8][1][10] = 3;
	auto test = arr;
	EXPECT_EQ(test[2][8][1][10], 3);
	EXPECT_EQ(test[9][7][8][12], 8);
}