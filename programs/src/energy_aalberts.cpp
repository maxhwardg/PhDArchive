//
// Created by max on 10/23/18.
//

#include "cxxopts.hpp"
#include "scorers/aalberts_scorer.hpp"

#include <string>
#include <iostream>

using namespace std;

int main(int argc, char **argv) {
    cxxopts::Options
            options("Energy Calculator (Aalberts & Nandagopal)",
                    "Calculates the free energy change of primary sequence and secondary structure pairs. "
                    "Expects a sequence of primary sequence and secondary structure strings on standard in. "
                    "Each primary sequence should be a string containing only 'A', 'U', 'G' and 'C' characters. "
                    "Each secondary structure should be in dot-bracket notation.");

    options.add_options()
            ("d,data_path", "Path to data_tables folder", cxxopts::value<string>()->default_value("data_tables/"))
            ("a,lengtha", "The distance for length-a segments", cxxopts::value<double>()->default_value("6.2"))
            ("b,lengthb", "The distance for length-b segments", cxxopts::value<double>()->default_value("15"))
            ("C,Cval", "The additive constant C applied to the cost of a multi-loop",
             cxxopts::value<librnary::kcalmol_t>()->default_value("0"))
            ("v,verbose", "Prints a full breakdown of the energy calculation (including accoutrements)")
            ("h,help", "Print help");

    string data_tables;
    double lengtha, lengthb;
    librnary::kcalmol_t C;
    bool verbose = false;

    try {
        options.parse(argc, argv);
        data_tables = options["data_path"].as<string>();
        lengtha = options["lengtha"].as<double>();
        lengthb = options["lengthb"].as<double>();
        C = options["Cval"].as<librnary::kcalmol_t>();
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

    librnary::AalbertsModel model(data_tables);
    model.SetNCoeffBase(lengtha);
    model.SetMCoeffBase(lengthb);
    model.SetAdditiveConstant(C);

    librnary::AalbertsScorer scorer(model);


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