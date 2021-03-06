//
// Created by max on 6/20/16.
//

#include "rna_structure.hpp"
#include "read_cts.hpp"

using namespace std;

librnary::CTData librnary::ReadCTFile(const string &full_path) {
	auto struc = librnary::LoadStructure(full_path);
	auto prim = librnary::StructureToPrimary(*struc);
	auto ss = librnary::StructureToMatching(*struc);
	return {full_path, prim, ss};
}

librnary::CTData librnary::ReadCTFile(const string &directory, const string &file) {
	auto struc = librnary::LoadStructure(directory + "/" + file);
	auto prim = librnary::StructureToPrimary(*struc);
	auto ss = librnary::StructureToMatching(*struc);
	return {file, prim, ss};
}


vector<librnary::CTData> librnary::ReadAllCTs(const string &base_path, istream &in) {
	auto sets = ReadFilesInCTSetFormat(base_path, in);
	assert(sets.size() == 1);
	vector<librnary::CTData> result;
	for (const auto& set : sets) {
		for (const auto& ct : set) {
			result.push_back(ct);
		}
	}
	return result;
}


vector<vector<librnary::CTData>> librnary::ReadFilesInCTSetFormat(const string &base_path, istream &in) {
	// A list of all the data sets.
	// Start with a single empty data set.
	vector<vector<librnary::CTData>> data_sets(1);

	while (in.good()) {
		in >> ws;
		string fname;
		getline(in, fname);

		if (fname == "end")
			break;

		// Start a new data set
		if (fname == "end set") {
			data_sets.emplace_back();
			continue;
		}

		data_sets.back().push_back(ReadCTFile(base_path, fname));
		in >> ws;
	}
	return data_sets;
}
