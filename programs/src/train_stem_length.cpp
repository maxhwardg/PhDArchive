//
// Created by max on 10/23/18.
//


#include <string>
#include <sstream>
#include <scorers/stem_length_scorer.hpp>
#include "folders/stem_length_folder.hpp"
#include "training/generic_ibf_trainer.hpp"
#include "cxxopts.hpp"

using namespace std;

const int STEM_LENGTH_COST_VECTOR_SIZE = 5;

class StemLengthParamSet {
protected:
    vector<librnary::energy_t> length_costs;
public:
    explicit StemLengthParamSet(vector<librnary::energy_t> costs) : length_costs(move(costs)) {
        assert(length_costs.size() == STEM_LENGTH_COST_VECTOR_SIZE);
    }
    string to_string() const {
        stringstream ss;
        ss << "Stem Length Costs: ";
        for (auto cost : length_costs) {
            ss << cost << ", ";
        }
        return ss.str();
    }
    bool operator==(const StemLengthParamSet &rhs) const {
        return length_costs == rhs.length_costs;
    }
    void LoadInto(librnary::StemLengthModel &model) const {
        model.SetLengthCosts(length_costs);
    }
    class SSInfo {
    protected:
        std::vector<int> stem_size_count;
    public:
        SSInfo() = default;
        explicit SSInfo(const librnary::Matching &match) {
            stem_size_count.assign(STEM_LENGTH_COST_VECTOR_SIZE, 0);
            auto stems = librnary::ExtractStems(match);
            for (const auto &stem : stems) {
                if (stem.num_pairs > static_cast<int>(stem_size_count.size()))
                    stem_size_count.back() += 1;
                else
                    stem_size_count[stem.num_pairs - 1]++;
            }
        }
        librnary::energy_t EnergyCost(const librnary::StemLengthModel &model) const {
            librnary::energy_t sum = 0;
            for (int i = 0; i < static_cast<int>(stem_size_count.size()); ++i) {
                sum += stem_size_count[i] * model.StemLengthCost(i + 1);
            }
            return sum;
        }
    };
};

int main(int argc, char **argv) {
    cxxopts::Options
            options("Train Stem Length Penalties",
                    "Trains the parameters of stem length penalties using IBF. "
                    "Expects a .ctset file as input on standard in. "
                    "For some .ctset files, see the data_set directory.");

    options.add_options()
            ("d,data_path", "Path to data_tables", cxxopts::value<string>()->default_value("data_tables/"))
            ("c,ct_path", "Path to the folder of CTs", cxxopts::value<string>()->default_value("data_set/ct_files/"))
            ("t,threads",
             "Number of threads to use",
             cxxopts::value<int>()->default_value(std::to_string(std::thread::hardware_concurrency())))
            ("l,disable_lonely_pairs", "Give lonely pairs a big energy penalty")
            ("h,help", "Print help");

    string data_tables, ct_path;
    size_t threads;
    bool no_lonely_pairs = false;

    try {
        options.parse(argc, argv);
        data_tables = options["data_path"].as<string>();
        ct_path = options["ct_path"].as<string>();
        threads = static_cast<size_t>(options["threads"].as<int>());
        if (options.count("disable_lonely_pairs") == 1) {
            no_lonely_pairs = true;
        }
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
    vector<StemLengthParamSet> params;

    librnary::energy_t a_min = -10, a_max = 30;

    if (no_lonely_pairs) a_min = a_max = 50;

    for (librnary::energy_t a = a_min; a <= a_max; ++a) {
        for (librnary::energy_t b = -10; b <= 30; ++b) {
            for (librnary::energy_t c = -10; c <= 30; ++c) {
                for (librnary::energy_t d = -10; d <= 30; ++d) {
                    for (librnary::energy_t e = -10; e <= 30; ++e) {
                        params.emplace_back(vector<librnary::energy_t>({e, a, b, c, d}));
                    }
                }
            }
        }
    }

    // Set up needed instances.
    librnary::StemLengthModel model(data_tables);  // Zero cost stems by default.
    librnary::StemLengthFolder folder(model);
    folder.SetMaxTwoLoop(30);
    folder.SetLonelyPairs(true);

    // Make the trainer and train!
    librnary::GenericIBFTrainer<StemLengthParamSet,
            librnary::StemLengthModel,
            librnary::StemLengthScorer,
            librnary::StemLengthFolder> trainer(model, cts, cout);
    trainer.SetNumStructureSeeds(5);
    trainer.SetThreads(threads);
    auto best_params = trainer.Train(params, params.front(), folder, 50);

    cout << "Best parameters: " << best_params.to_string() << endl;


    return 0;
}