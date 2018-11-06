//
// Created by max on 8/16/17.
//

#include <string>
#include <sstream>

#include "cxxopts.hpp"

#include "read_cts.hpp"
#include "training/IBF_multiloop.hpp"
#include "models/asymmetry_model.hpp"
#include "folders/asymmetry_folder.hpp"
#include "scorers/asymmetry_scorer.hpp"


using namespace std;

class AsymmetryParamSet {
protected:
    librnary::energy_t init, branch, unpaired, asym_cost;
public:
    AsymmetryParamSet(librnary::energy_t ini,
                      librnary::energy_t br,
                      librnary::energy_t up, librnary::energy_t a_cost)
            : init(ini), branch(br), unpaired(up), asym_cost(a_cost) {}

    string to_string() const {
        stringstream ss;
        ss << "init = " << init << " branch = " << branch << " unpaired = " << unpaired << " asymmetry cost = "
           << asym_cost;
        return ss.str();
    }

    bool operator==(const AsymmetryParamSet &rhs) const {
        return init == rhs.init && branch == rhs.branch && unpaired == rhs.unpaired
               && asym_cost == rhs.asym_cost;
    }

    void LoadInto(librnary::AsymmetryModel &model) const {
        model.SetMLParams(init, branch, unpaired, asym_cost);
    }
    class MultiInfo {
    protected:
        int branches{}, unpaired{}, sum_asymmetry{};
    public:
        MultiInfo() = default;
        MultiInfo(const librnary::Surface &surf) {
            auto lr = librnary::LoopRegion(surf);
            branches = librnary::ExtractBranches(lr);
            unpaired = librnary::ExtractUnpaired(lr);
            sum_asymmetry = librnary::ExtractSumAsymmetry(lr);
        }
        librnary::energy_t MLClosure(const librnary::AsymmetryModel &model) const {
            return model.MLClosure(branches, unpaired, sum_asymmetry);
        }
    };
};

int main(int argc, char **argv) {
    cxxopts::Options
            options("Trains Linear Asymmetry Model",
                    "Trains the parameters of the linear asymmetry model using IBF. "
                    "Expects a .ctset file as input on standard in. "
                    "For some .ctset files, see the data_set directory.");

    options.add_options()
            ("d,data_path", "Path to data_tables", cxxopts::value<string>()->default_value("data_tables/"))
            ("c,ct_path", "Path to the folder of CTs", cxxopts::value<string>()->default_value("data_set/ct_files/"))
            ("t,threads",
             "Number of threads to use",
             cxxopts::value<int>()->default_value(std::to_string(std::thread::hardware_concurrency())))
            ("h,help", "Print help");

    string data_tables, ct_path;
    size_t threads;

    try {
        options.parse(argc, argv);
        data_tables = options["data_path"].as<string>();
        ct_path = options["ct_path"].as<string>();
        threads = static_cast<size_t>(options["threads"].as<int>());
        if (options.count("help") == 1) {
            cout << options.help({"", "Group"}) << endl;
            return 0;
        }

    } catch (const cxxopts::OptionException &e) {
        cout << "Argument parsing error: " << e.what() << endl;
        return 1;
    }

    // Read CTs.

    auto cts = librnary::ReadFilesInCTSetFormat(ct_path, cin);

    // Generate parameter list.
    vector<AsymmetryParamSet> params;

    for (librnary::energy_t init = 40; init <= 205; ++init) {
        for (librnary::energy_t br = -35; br <= 30; ++br) {
            for (librnary::energy_t up = -35; up <= 30; ++up) {
                for (librnary::energy_t asym_cost = -30; asym_cost <= 30; ++asym_cost) {
                    params.emplace_back(init, br, up, asym_cost);
                }
            }
        }
    }

    // Set up needed instances.
    librnary::AsymmetryModel model(data_tables);
    librnary::AsymmetryFolder folder(model);
    folder.SetMaxTwoLoop(30);
    folder.SetLonelyPairs(false);
    // Give multi-loops no cost.
    model.SetMLParams(0, 0, 0, 0);

    // Make the trainer and train!
    librnary::IBFMultiLoop<AsymmetryParamSet,
            librnary::AsymmetryModel,
            librnary::AsymmetryScorer,
            librnary::AsymmetryFolder> trainer(model, cts, clog);
    trainer.SetNumStructureSeeds(5);
    trainer.SetThreads(threads);
    auto best_params = trainer.Train(params, params.front(), folder, 50);

    cout << "Best parameters: " << best_params.to_string() << endl;


    return 0;
}