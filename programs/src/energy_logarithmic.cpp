//
// Created by max on 10/23/18.
//

#include "cxxopts.hpp"
#include "scorers/nn_scorer.hpp"
#include "models/nn_unpaired_model.hpp"

#include <string>
#include <iostream>

using namespace std;

int main(int argc, char **argv) {
    cxxopts::Options
            options("Energy Calculator (logarithmic)",
                    "Calculates the free energy change of primary sequence and secondary structure pairs. "
                    "Expects a sequence of primary sequence and secondary structure strings on standard in. "
                    "Each primary sequence should be a string containing only 'A', 'U', 'G' and 'C' characters. "
                    "Each secondary structure should be in dot-bracket notation.");

    options.add_options()
            ("d,data_path", "Path to data_tables folder", cxxopts::value<string>()->default_value("data_tables/"))
            ("i,ml_init",
             "Initiation cost (in tenth kcal/mol) of multi-loop (A in A+B*branches+C*unpaired_up_to_E+*D*ln(unpaired/E))",
             cxxopts::value<librnary::energy_t>()->default_value("101"))
            ("b,ml_branch",
             "Branch cost (in tenth kcal/mol) in a multi-loop (B in A+B*branches+C*unpaired_up_to_E+*D*ln(unpaired/E))",
             cxxopts::value<librnary::energy_t>()->default_value("-3"))
            ("u,ml_unpaired",
             "Unpaired cost (in tenth kcal/mol) in a multi-loop (C in A+B*branches+C*unpaired_up_to_E+*D*ln(unpaired/E))",
             cxxopts::value<librnary::energy_t>()->default_value("-3"))
            ("m,ml_log_mult", "Logarithmic coefficient (in kcal/mol for extra accuracy) in multi-loop model "
                              "(D in A+B*branches+C*unpaired_up_to_E+*D*ln(unpaired/E))",
             cxxopts::value<librnary::kcalmol_t>()->default_value("1.1"))
            ("p,ml_pivot", "Pivot (in number of unpaired) to a logarithmic term in multi-loop model "
                           "(E in A+B*branches+C*unpaired_up_to_E+*D*ln(unpaired/E))",
             cxxopts::value<librnary::energy_t>()->default_value("6"))
            ("v,verbose", "Prints a full breakdown of the energy calculation (including accoutrements)")
            ("h,help", "Print help");

    string data_tables;
    librnary::energy_t ml_init, ml_branch, ml_unpaired;
    librnary::kcalmol_t ml_log_mult;
    int ml_pivot;
    bool verbose = false;

    try {
        options.parse(argc, argv);
        data_tables = options["data_path"].as<string>();
        ml_init = options["ml_init"].as<librnary::energy_t>();
        ml_branch = options["ml_branch"].as<librnary::energy_t>();
        ml_unpaired = options["ml_unpaired"].as<librnary::energy_t>();
        ml_log_mult = options["ml_log_mult"].as<librnary::kcalmol_t>();
        ml_pivot = options["ml_pivot"].as<int>();
        if (options.count("verbose") == 1) {
            verbose = true;
        }
        if (options.count("help") == 1) {
            cout << options.help({"", "Group"}) << endl;
            return 0;
        }

    } catch (const cxxopts::OptionException &e) {
        cout << "Argument parsing error: " << e.what() << endl;
        return 1;
    }

    librnary::NNUnpairedModel model(data_tables);
    model.SetMLInitConstant(ml_init);
    model.SetMLBranchCost(ml_branch);
    model.SetMLUnpairedCost(ml_unpaired);
    model.SetMLLogMultiplier(ml_log_mult);
    model.SetMLUnpairedPivot(ml_pivot);

    librnary::NNScorer<librnary::NNUnpairedModel> scorer(model);


    string primary_str, dot_bracket;
    while (cin >> primary_str >> dot_bracket) {
        auto primary = librnary::StringToPrimary(primary_str);
        scorer.SetRNA(primary);
        librnary::SSTree ss_tree(librnary::DotBracketToMatching(dot_bracket));
        librnary::energy_t e = scorer.ScoreExterior(ss_tree.RootSurface());
        cout << "Free Energy Change: " << librnary::EnergyToKCal(e) << " (kcal/mol)" << endl;
        if (verbose) {
            auto surface_score = scorer.TraceExterior(ss_tree.RootSurface());
            cout << surface_score.Describe(' ', 0) << endl;
        }
    }
}