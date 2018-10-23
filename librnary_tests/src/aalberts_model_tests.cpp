//
// Created by max on 6/22/16.
//

#include <gtest/gtest.h>
#include "models/aalberts_model.hpp"

using namespace std;

const string DATA_TABLE_PATH = "../../data_tables/";

TEST(AalbertsModel, Eqn3FromPaper) {
	librnary::AalbertsModel model(DATA_TABLE_PATH);
	// Gfjc2(11,0) - Gfjc2(12,4) == -1
	EXPECT_EQ(model.MLInit(11, 0) - model.MLInit(12, 4), -10);
	EXPECT_EQ(model.MLInitUpBr(11, 0) - model.MLInitUpBr(8, 4), -10);
}