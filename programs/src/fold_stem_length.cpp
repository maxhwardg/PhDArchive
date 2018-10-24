//
// Created by max on 10/23/18.
//

#include "cxxopts.hpp"
#include "folders/stem_length_folder.hpp"

#include <string>
#include <iostream>

using namespace std;

int main(int argc, char **argv) {
    cxxopts::Options
            options("RNA Folding Algorithm (stem length penalties)",
                    "Takes an RNA primary sequence and applies a folding algorithm to predict a MFE structure. "
                    "Uses the classic linear (or affine) model of multi-loop free energy. "
                    "Include penalties for stem length. "
                    "Reads a list of space separated list of RNA primary sequences on standard input"
                    " and folds each in turn."
                    "RNA sequence should be a string of only 'A', 'G', 'U', or 'C' characters."
                    "Outputs the results to standard out.");

    options.add_options()
            ("d,data_path", "Path to data_tables folder", cxxopts::value<string>()->default_value("data_tables/"))
            ("i,ml_init", "Initiation cost (in tenth kcal/mol) of multi-loop (A in A+B*branches+C*unpaired)",
             cxxopts::value<librnary::energy_t>()->default_value("93"))
            ("b,ml_branch", "Branch cost (in tenth kcal/mol) in a multi-loop (B in A+B*branches+C*unpaired)",
             cxxopts::value<librnary::energy_t>()->default_value("-6"))
            ("u,ml_unpaired", "Unpaired cost (in tenth kcal/mol) in a multi-loop (C in A+B*branches+C*unpaired)",
             cxxopts::value<librnary::energy_t>()->default_value("0"))
            ("s,stem_length_costs", "A string that is a space separated cost of stems of increasing lengths. "
                                    "For example, '50 6 15 15 9' "
                                    "to have stems of length 1 cost 5.0, length 2 cost 0.6, and so on.",
             cxxopts::value<string>()->default_value("50 6 15 15 9"))
            ("t,two_loop_max_size",
             "The maximum number of unpaired nucleotides allowed in a two-loop. "
             "(Set this to a very large number for unlimited)",
             cxxopts::value<int>()->default_value("30"))
            ("h,help", "Print help");

    string data_tables;
    librnary::energy_t ml_init, ml_branch, ml_unpaired;
    int max_two_loop_size;
    vector<librnary::energy_t> stem_length_costs;

    try {
        options.parse(argc, argv);
        data_tables = options["data_path"].as<string>();
        ml_init = options["ml_init"].as<librnary::energy_t>();
        ml_branch = options["ml_branch"].as<librnary::energy_t>();
        ml_unpaired = options["ml_unpaired"].as<librnary::energy_t>();
        max_two_loop_size = options["two_loop_max_size"].as<int>();

        std::stringstream ss(options["stem_length_costs"].as<string>());
        librnary::energy_t e;
        while (ss >> e) {
            stem_length_costs.push_back(e);
        }

        if (options.count("help") == 1) {
            cout << options.help({"", "Group"}) << endl;
            return 0;
        }

    } catch (const cxxopts::OptionException &e) {
        cout << "Argument parsing error: " << e.what() << endl;
        return 1;
    }

    librnary::StemLengthModel model(data_tables);
    model.SetMLInitCost(ml_init);
    model.SetMLBranchCost(ml_branch);
    model.SetMLUnpairedCost(ml_unpaired);
    model.SetLengthCosts(stem_length_costs);

    librnary::StemLengthFolder folder(model);
    folder.SetMaxTwoLoop(max_two_loop_size);

    string primary_str;
    while (cin >> primary_str) {
        auto primary = librnary::StringToPrimary(primary_str);
        librnary::energy_t e = folder.Fold(primary);
        cout << librnary::MatchingToDotBracket(folder.Traceback()) << endl;
        cout << "MFE: " << librnary::EnergyToKCal(e) << " (kcal/mol)" << endl;
    }
}