//
// Created by max on 8/29/16.
//

#include <gtest/gtest.h>

#include <energy.hpp>

TEST(ENERGY, ConvertKcalToEnergy) {
	EXPECT_EQ(1, librnary::KCalToEnergy(0.1));
	EXPECT_EQ(-45, librnary::KCalToEnergy(-4.53));
	EXPECT_EQ(114, librnary::KCalToEnergy(11.395));
}

TEST(ENERGY, ConvertEnergyToKcal) {
	EXPECT_EQ(0.1, librnary::EnergyToKCal(1));
	EXPECT_EQ(-4.5, librnary::EnergyToKCal(-45));
	EXPECT_EQ(11.4, librnary::EnergyToKCal(114));
}

// This is set this way to provide information about the above tests failing due to a changed Iota.
TEST(ENERGY, ConverstionFactor) {
	EXPECT_EQ(librnary::EnergyIota, 0.1);
}