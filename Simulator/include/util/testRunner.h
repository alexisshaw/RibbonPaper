//
// Created by 13114347 on 11/10/2020.
//

#ifndef STEANE_SIM_TESTRUNNER_H
#define STEANE_SIM_TESTRUNNER_H

#include <iostream>
#include <chrono>
#include <vector>
#include <iterator>
#include "sim/chp.h"
#include "codes/Steane.h"
#include "sim/balanced_error_sim.h"
#include "stats.hpp"
#include "prettyprint.h"

using std::chrono::high_resolution_clock;
using std::chrono::microseconds;
using std::chrono::duration_cast;

using CodeTools::BalancedErrorSim;
using CodeTools::ChpSimulator;
using CodeTools::single_pauli_state;

template<typename SimBase, template<typename baseSim> typename SimT, typename Functor>
class testRunner {
    Functor f;
public:
    explicit testRunner(Functor f_in): f{f_in}{}

    auto getSim(long double p){
        std::shared_ptr<std::mt19937_64> rand = getRandGen();
        ChpSimulator<SimBase::REQUIRED_QUBITS + SimBase::REQUIRED_ANCS> chp = ChpSimulator<SimBase::REQUIRED_QUBITS + SimBase::REQUIRED_ANCS>(rand);
        BalancedErrorSim<decltype(chp)> balSim = BalancedErrorSim(p,p,0,0,chp, rand);
        return balSim;
    }

    void runSimulations(int n, int print_interval, const std::vector<long double> &error_rates) {
        auto start = high_resolution_clock::now();
        auto current = high_resolution_clock::now();
        for(auto p: error_rates){
            double duration = duration_cast<microseconds>(current - start).count() / 1E6;
            std::cout << duration << " ";

            GetRunStats(p, n, print_interval);

            current = high_resolution_clock::now();
        }

        double duration = duration_cast<microseconds>(current - start).count() / 1E6;
        std::cout << duration << " " << std::endl;
    }

    std::vector<long double> getErrorSamples(int spacing, int decades, int startno, long double initial_rate) {
        std::vector<long double> error_rates{};
        for(int i=startno; i< (spacing*decades + 1) ; i++){
            long double il = i;
            long double spacingl = spacing;
            error_rates.push_back(powl(10.0L, -(il / spacingl))* initial_rate);
        }
        return error_rates;
    }

    void GetRunStats(long double p, int n, int print_interval) {
        std::vector<uint_fast64_t> results;
        auto start = high_resolution_clock::now();

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
            uint_fast64_t len = get_fault_error_syndromes(sim, SimBase::getProductInitialState(single_pauli_state(i % 6)), max);
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

        #pragma omp master
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
    uint_fast64_t get_fault_error_syndromes(const ChpSim &sim, const typename SimBase::InitialState& state, uint_fast64_t max_iterations){
        SimT<ChpSim> code_sim{sim, state};
        uint_fast64_t retval = max_iterations;
        for (uint_fast64_t i = 0; i < max_iterations; i++){
            if(f(code_sim)){
                retval = i;
                break;
            }
        }
        return retval;
    }
};


#endif //STEANE_SIM_TESTRUNNER_H
