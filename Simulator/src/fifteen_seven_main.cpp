#include "util/FTCorrectorFactory.h"
#include "codes/FifteenSeven.h"
#include "util/testRunner.h"
#include "absl/container/flat_hash_map.h"

using CodeTools::FifteenSeven::FifteenSevenBase;
using CodeTools::FifteenSeven::FifteenSeven;

int main() {
    auto factory = CorrectorFactory::FTCorrectorFactory<FifteenSevenBase, FifteenSeven>();
    auto map =  factory.getAndTestCorrectionMap();
    CodeTools::FifteenSeven::FifteenSevenBase::correction_map = map;

    int n = 900;
    int print_interval = 1000;

    int spacing = 6;
    int decades = 7;
    int startno = 1;
    long double initial_rate = 1E-1L;

    auto error_correction_round = [](auto &fifteenSim) -> bool{
        fifteenSim.correctFTSyndrome(fifteenSim.extractFTSyndrome());
        return !fifteenSim.confirm_no_error();
    };

    testRunner<FifteenSevenBase, FifteenSeven, decltype(error_correction_round)> runner(error_correction_round);

    std::vector<long double> error_rates = runner.getErrorSamples(spacing, decades, startno,initial_rate);
    std::cout << "p: " << error_rates << std::endl;

    runner.runSimulations(n, print_interval, error_rates);

    return 0;
}
