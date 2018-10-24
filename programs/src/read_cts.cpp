//
// Created by max on 10/23/18.
//

#include "cxxopts.hpp"
#include "read_cts.hpp"

#include <string>
#include <iostream>

using namespace std;

int main(int argc, char **argv) {
    cxxopts::Options
            options("CT Reader/Parser",
                    "Takes a sequence of space separated CT paths on standard in and produces extracted primary "
                    "sequences and secondary structures on standard out.\n\n"
                    "Useful to pipe the result of this program to a fold/efn program.");

    options.add_options()
            ("s,secondary_structure", "Setting this flag will cause the program to output the secondary structure")
            ("h,help", "Print help");

    bool secondary_structure = false;

    try {
        options.parse(argc, argv);
        if (options.count("secondary_structure") == 1) {
            secondary_structure = true;
        }
        if (options.count("help") == 1) {
            cout << options.help({"", "Group"}) << endl;
            return 0;
        }

    } catch (const cxxopts::OptionException &e) {
        cout << "Argument parsing error: " << e.what() << endl;
        return 1;
    }

    string ct_path;
    while (cin >> ct_path) {
        auto ct = librnary::ReadCTFile(ct_path);
        auto primary = librnary::PrimaryToString(ct.primary);
        cout << primary << endl;
        if (secondary_structure) {
            auto secondary = librnary::MatchingToDotBracket(ct.match);
            cout << secondary << endl;
        }
    }
}