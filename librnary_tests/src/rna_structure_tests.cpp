//
// Created by max on 6/20/16.
//

/*
 * Contains tests for the rna_structure wrapper header.
 * These tests basically ensure internal consistency and do some sanity checks.
 */

#include <gtest/gtest.h>
#include <scorers/nn_scorer.hpp>
#include <models/nn_model.hpp>
#include "models/nn_affine_model.hpp"
#include <ss_enumeration.hpp>

#include "random.hpp"

using namespace std;

const string DATA_TABLE_PATH = "../../data_tables/";

TEST(RNAstructure, FuzzConversionToStructureAndBack) {
	auto re = librnary::RandomEngineForTests();
	const unsigned CASES = 500, RNA_LEN = 50, SS_TRIALS = 50;
	for (unsigned tc = 0; tc < CASES; ++tc) {
		auto prim = librnary::RandomPrimary(re, RNA_LEN);
		auto match = librnary::RandomMatching(prim, re, SS_TRIALS);
		auto struc = librnary::LoadStructure(prim, match);
		EXPECT_EQ(prim, librnary::StructureToPrimary(*struc));
		EXPECT_EQ(match, librnary::StructureToMatching(*struc));
	}
}

TEST(RNAstructure, FoldAUCGUAGCUUAUUAGCGGCGAGCUACGUACGGUACGUAC) {
	auto prim = librnary::StringToPrimary("AUCGUAGCUUAUUAGCGGCGAGCUACGUACGGUACGUAC");
	auto struc = librnary::LoadStructure(prim);
	auto dt = librnary::LoadDatatable(DATA_TABLE_PATH);
	EXPECT_EQ(librnary::RunMFEFold(*dt, *struc), -114);
	EXPECT_EQ(librnary::RunEFN2(*dt, *struc), -114);
	EXPECT_EQ(librnary::RunEFN2WithSimpleMulti(*dt, *struc), -114);
}

TEST(RNAstructure, FoldLargeMultiLoop) {
	auto prim = librnary::StringToPrimary(
		"UCGAUCUACGUGUCGUGAUCGUCGUCGUAGCUAGCUGACUGAUGAUGUAGCUAGCUAGCUAGCUAGUGAUCGAUGUACGUAC"
			"GUCAGAUCGAUCGCUGACUGAUAGUACGUCGAUUACGUACGUAGCUGCUGACUAGCUGCAUCGUAGCUGAUGCUAGUCGAUCGUACGUCAGUCGAUGC"
			"UAGCUAGCGCUAUAUUAUCUGACGAUCUAUCGAGCGACGGCGCUAUGCUAGUCGUGACGUAUCGUAUGCUGAUCGAGGCGGCGAGCGAUUA"
			"UAUAUAUAGCGCGAUCUUACUAUCGACUGUAGCGAUCGUAGCUAUCUGACGUACGUACGACUGUACAUCGUAGUACGUACG");
	auto struc = librnary::LoadStructure(prim);
	auto dt = librnary::LoadDatatable(DATA_TABLE_PATH);
	EXPECT_EQ(librnary::RunMFEFold(*dt, *struc), -1439);
	EXPECT_EQ(librnary::RunEFN2(*dt, *struc), -1427);
	EXPECT_EQ(librnary::RunEFN2WithSimpleMulti(*dt, *struc), -1439);
}