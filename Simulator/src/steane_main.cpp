#include <iostream>
#include <vector>
#include "util/testRunner.h"

using CodeTools::Steane::SteaneSimBase;
using CodeTools::Steane::SteaneSimT;

int main() {
    int n = 900;
    int print_interval = 1000;

    int spacing = 12;
    int decades = 5;
    int startno = 1;
    long double initial_rate = 1E-1L;

    auto error_correction_round = [](auto &steaneSim) -> bool{
        steaneSim.correctFTSyndrome(steaneSim.extractFTSyndrome());
        return !steaneSim.confirm_no_error();
    };

    testRunner<CodeTools::Steane::SteaneSimBase, CodeTools::Steane::SteaneSimT, decltype(error_correction_round)> runner(error_correction_round);

    std::vector<long double> error_rates = runner.getErrorSamples(spacing, decades, startno,initial_rate);
    std::cout << "p: " << error_rates << std::endl;

    runner.runSimulations(n, print_interval, error_rates);

    return 0;
}

