#include <iostream>
#include <chrono>
#include <vector>
#include <iterator>
#include "chp.h"
#include "SteaneSimBase.h"
#include "balanced_error_sim.h"
#include "stats.hpp"

using namespace std::chrono;
using namespace SteaneSim;
using namespace SteaneSim::Steane;

template<typename ChpSim>
uint_fast64_t get_fault_error_syndromes(const ChpSim &sim, const InitialState& state, uint_fast64_t max_iterations);
std::shared_ptr<std::mt19937_64> getRandGen();

auto getSim(long double p){
    std::shared_ptr<std::mt19937_64> rand = getRandGen();
    ChpSimulator<REQUIRED_QUBITS + REQUIRED_ANCS> chp = ChpSimulator<REQUIRED_QUBITS + REQUIRED_ANCS>(rand);
    BalancedErrorSim<decltype(chp)> balSim = BalancedErrorSim(p,p,0,0,chp, rand);
    return balSim;
}

void GetRunStats(long double p, int n, int print_interval);

template <typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
    if ( !v.empty() ) {
        out << '[';
        std::copy (v.begin(), v.end(), std::ostream_iterator<T>(out, ", "));
        out << "\b\b]";
    }
    return out;
}

int main() {
    int n = 900;
    int print_interval = 1000;

    std::vector<long double> error_rates{};
    int spacing = 12;
    int decades = 5;
    int startno = 1;
    for(int i=startno; i< (spacing*decades + 1) ; i++){
        long double il = i;
        long double spacingl = spacing;
        error_rates.push_back(powl(10.0L, -(il / spacingl))* 1E-1L);
    }
    std::cout << "p: " << error_rates << std::endl;
    auto start = high_resolution_clock::now();
    for(auto p: error_rates){
        auto current = high_resolution_clock::now();
        double duration = duration_cast<microseconds>(current - start).count() / 1E6;
        std::cout << duration << " ";
        GetRunStats(p, n, print_interval);
    }
    auto current = high_resolution_clock::now();
    double duration = duration_cast<microseconds>(current - start).count() / 1E6;
    std::cout << duration << " " << std::endl;
    return 0;
}

void GetRunStats(long double p, int n, int print_interval) {
    std::vector<uint_fast64_t> results;

    // This...
    #pragma omp master
    {
        results = std::vector<uint_fast64_t>{};
    }
    #pragma omp barrier

    auto total_start = high_resolution_clock::now();

    #pragma omp parallel for schedule(dynamic) shared(results, std::cout) firstprivate(total_start, p, n, print_interval) default(none)
    for(int i = 0; i < n; ++i) {
        auto start = high_resolution_clock::now();
        auto sim = getSim(p);
        uint_fast64_t max = uint_fast64_t(1E11);
        uint_fast64_t len = get_fault_error_syndromes(sim, single_pauli_state(i % 6), max);
        #pragma omp critical(results)
        {
            results.push_back(len);
            auto end = high_resolution_clock::now();
            if(len == max) std::cout << "Max reached for index: "<< i << std::endl;
            if((i%print_interval) == print_interval-1){
                std::cout << "iteration no: " << i
                          << " time: "<< duration_cast<microseconds>(end - start).count()/1E6
                          << " elapsed: "<< duration_cast<microseconds>(end - total_start).count()/1E6
                          << " len: " << len << std::endl;
        };
        }
    }

    #pragma omp single
    {
        auto total_end = high_resolution_clock::now();
        double duration = duration_cast<microseconds>(total_end - total_start).count() / 1E6;
        double mean_value = mean(results);
        double std_dev = standarddeviation(results);
        double median_value = median(results);
        std::cout << "Probability " << p << " stats:" << "  num runs: " << n << "  time: " << duration
                  << "  Mean: " << mean_value << " ~ " << 1 / mean_value << " logical" << "  std deviation: "
                  << std_dev  << "  Median: " << median_value << std::endl;
    }
}

std::shared_ptr<std::mt19937_64> getRandGen(){
    std::random_device rd{"/dev/random"};
    std::seed_seq seq{rd(), rd(), rd(), rd()};
    std::shared_ptr<std::mt19937_64> gen = std::make_shared<std::mt19937_64>(seq);
    return gen;
}

template<typename ChpSim>
uint_fast64_t get_fault_error_syndromes(const ChpSim &sim, const InitialState& state, uint_fast64_t max_iterations){
    SteaneSimBase<ChpSim> steane_sim{sim, state};
    for (uint_fast64_t i = 0; i < max_iterations; i++){
        auto a1 = steane_sim.perform_round();
        auto a2 = steane_sim.perform_round();
        if (a1 == a2){
            steane_sim.simple_correct(a1);
        } else {
            auto a3 = steane_sim.perform_round();
            if (a2 == a3){
                steane_sim.fault_tolerance_correct(a1, a2);
            } else {
                steane_sim.fault_tolerance_correct(a2, a3);
            }
        }
        if(!steane_sim.confirm_no_error())
            return i;
    }
    return max_iterations;
}
