//
// Created by max on 6/2/16.
//

#include <gtest/gtest.h>
#include "models/nn_unpaired_model.hpp"
#include "models/nn_affine_model.hpp"

using namespace std;

const string DATA_TABLE_PATH = "../../data_tables/";

// All these tests use NNAffineModel because it is the simplest model that implements the virtual NNModel class.

TEST(NNModel, Hairpin_AAAAAAAU) {
	auto primary = librnary::StringToPrimary("AAAAAAAU");
	librnary::NNAffineModel nn_model(DATA_TABLE_PATH, primary);
	EXPECT_EQ(nn_model.OneLoop(0, primary.size() - 1), 54 + 5 - 8);
}

TEST(NNModel, Hairpin_AGGAAGU) {
	auto primary = librnary::StringToPrimary("AGGAAGU");
	librnary::NNAffineModel nn_model(DATA_TABLE_PATH, primary);
	EXPECT_EQ(nn_model.OneLoop(0, primary.size() - 1), 57 - 8 - 8 + 5);
}

TEST(NNModel, Hairpin_CCGAGG) {
	auto primary = librnary::StringToPrimary("CCGAGG");
	librnary::NNAffineModel nn_model(DATA_TABLE_PATH, primary);
	EXPECT_EQ(nn_model.OneLoop(0, primary.size() - 1), 35);
}

TEST(NNModel, Hairpin_ACCCCCCU) {
	auto primary = librnary::StringToPrimary("ACCCCCCU");
	librnary::NNAffineModel nn_model(DATA_TABLE_PATH, primary);
	EXPECT_EQ(nn_model.OneLoop(0, primary.size() - 1), 16 + 6 * 3 + 54 - 7 + 5);
}

TEST(NNModel, Hairpin_CGGGGGAAGUCCG) {
	auto primary = librnary::StringToPrimary("CGGGGGAAGUCCG");
	librnary::NNAffineModel nn_model(DATA_TABLE_PATH, primary);
	EXPECT_EQ(nn_model.OneLoop(3, 9), -22 + 57 - 8 - 8 + 5);
}

TEST(NNModel, Stack_AGCGCUAGCGCU) {
	auto primary = librnary::StringToPrimary("AGCGCUAGCGCU");
	librnary::NNAffineModel nn_model(DATA_TABLE_PATH, primary);
	EXPECT_EQ(nn_model.TwoLoop(0, 1, primary.size() - 2, primary.size() - 1), -21);
	EXPECT_EQ(nn_model.TwoLoop(1, 2, primary.size() - 3, primary.size() - 2), -34);
	EXPECT_EQ(nn_model.TwoLoop(2, 3, primary.size() - 4, primary.size() - 3), -24);
}

TEST(NNModel, Bulge_GAACAGAAACUC) {
	auto primary = librnary::StringToPrimary("GAACAGAAACUC");
	//									      .(...(...)).
	librnary::NNAffineModel nn_model(DATA_TABLE_PATH, primary);
	EXPECT_EQ(nn_model.TwoLoop(1, 5, 9, 10), 5 + 32);
}

TEST(NNModel, Bulge_SingleC) {
	auto primary = librnary::StringToPrimary("GCCCGAAACGGC");
	//									      (((.(...))))
	librnary::NNAffineModel nn_model(DATA_TABLE_PATH, primary);
	// In RNAstructure version 5.8.1, this is off by 1 since the software truncates instead of rounding the bonus.
	EXPECT_EQ(nn_model.TwoLoop(2, 4, 8, 9), 38 - 9 - 24 - 7);
}

TEST(NNModel, Interal_2x2) {
	auto primary = librnary::StringToPrimary("CAGACGCGGAUA");
	//									      .(..(..)..).
	librnary::NNAffineModel nn_model(DATA_TABLE_PATH, primary);
	EXPECT_EQ(nn_model.TwoLoop(1, 4, 7, 10), -11);
}

TEST(NNModel, Interal_1x5) {
	auto primary = librnary::StringToPrimary("CAGCGCGGAAAGUG");
	//									      .(.(..).....).
	librnary::NNAffineModel nn_model(DATA_TABLE_PATH, primary);
	EXPECT_EQ(nn_model.TwoLoop(1, 3, 6, 12), 20 + 24 + 7);
}

TEST(NNModel, Interal_2x3) {
	auto primary = librnary::StringToPrimary("CAGACGCGGAGUG");
	//									      .(..(..)...).
	librnary::NNAffineModel nn_model(DATA_TABLE_PATH, primary);
	// This +1 is there because RNAstructure 5.8.1 has an off by 1 in its terminal mismatch parameters for 2x3 ILs.
	EXPECT_EQ(nn_model.TwoLoop(1, 4, 7, 11), (20 + 6 - 8 - 12 + 7) + 1);
}


TEST(NNModel, FlushCoax_AUCGG) {
	auto primary = librnary::StringToPrimary("AUCGG");
	librnary::NNAffineModel nn_model(DATA_TABLE_PATH, primary);
	EXPECT_EQ(nn_model.FlushCoax(0, 1, 2, 3), -24);
}

TEST(NNModel, MismatchCoax_AGCGGU) {
	auto primary = librnary::StringToPrimary("AGCGGU");
	librnary::NNAffineModel nn_model(DATA_TABLE_PATH, primary);
	EXPECT_EQ(nn_model.MismatchCoax(0, 5, 2, 3), -30);
}