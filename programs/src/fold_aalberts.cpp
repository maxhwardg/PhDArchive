//
// Created by max on 10/22/18.
//


#include <iostream>

#include "cxxopts.hpp"
#include "folders/aalberts_folder.hpp"

using namespace std;

int main(int argc, char **argv) {
    cxxopts::Options
            options("RNA Folding Algorithm (Aalberts & Nandagopal)",
                    "Takes an RNA primary sequence and applies a folding algorithm to predict a MFE structure. "
                    "Uses the Aalberts & Nandagopal model of multi-loop free energy. "
                    "Reads a list of space separated list of RNA primary sequences on standard input"
                    " and folds each in turn."
                    "RNA sequence should be a string of only 'A', 'G', 'U', or 'C' characters."
                    "Outputs the results to standard out.");

    options.add_options()
            ("d,data_path", "Path to data_tables folder", cxxopts::value<string>()->default_value("data_tables/"))
            ("a,lengtha", "The distance for length-a segments", cxxopts::value<double>()->default_value("6.2"))
            ("b,lengthb", "The distance for length-b segments", cxxopts::value<double>()->default_value("15"))
            ("C,Cval", "The additive constant C applied to the cost of a multi-loop",
             cxxopts::value<librnary::kcalmol_t>()->default_value("0"))
            ("t,two_loop_max_size",
             "The maximum number of unpaired nucleotides allowed in a two-loop. "
             "(Set this to a very large number for unlimited)",
             cxxopts::value<int>()->default_value("30"))
            ("l,lonely_pairs", "Setting this flag will disable the no lonely pairs heuristic")
            ("h,help", "Print help");

    string data_tables;
    double lengtha, lengthb;
    librnary::kcalmol_t C;
    bool lonely_pairs = false;
    int max_two_loop_size;

    try {
        options.parse(argc, argv);
        data_tables = options["data_path"].as<string>();
        lengtha = options["lengtha"].as<double>();
        lengthb = options["lengthb"].as<double>();
        C = options["Cval"].as<librnary::kcalmol_t>();
        max_two_loop_size = options["two_loop_max_size"].as<int>();
        if (options.count("lonely_pairs") == 1) {
            lonely_pairs = true;
        }
        if (options.count("help") == 1) {
            cout << options.help({"", "Group"}) << endl;
            return 0;
        }

    } catch (const cxxopts::OptionException &e) {
        cout << "Argument parsing error: " << e.what() << endl;
        return 1;
    }

    librnary::AalbertsModel model(data_tables);
    model.SetNCoeffBase(lengtha);
    model.SetMCoeffBase(lengthb);
    model.SetAdditiveConstant(C);

    librnary::AalbertsFolder folder(model);
    folder.SetMaxTwoLoop(max_two_loop_size);
    folder.SetLonelyPairs(lonely_pairs);

    string primary_str;
    while (cin >> primary_str) {
        auto primary = librnary::StringToPrimary(primary_str);
        librnary::energy_t e = folder.Fold(primary);
        cout << librnary::MatchingToDotBracket(folder.Traceback()) << endl;
        cout << "MFE: " << librnary::EnergyToKCal(e) << " (kcal/mol)" << endl;
    }
}