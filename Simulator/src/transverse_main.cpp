//
//  main.c
//  QuBOX
//
//  Created by Simon Devitt on 2/8/19.
//  Copyright © 2019 Simon Devitt. All rights reserved.
//

#include <iostream>
#include <sstream>
#include <random>
#include <boost/container/small_vector.hpp>
#include <boost/program_options.hpp>
#include <utility>
#include <numbers>
#include <numeric>
#include <omp.h>
#include <chrono>

#include "sim/state_vector.h"
#include "sim/auto_balanced_error_sim.h"
#include "sim/chp.h"
#include "codes/Planar.h"
#include "prettyprint.h"
#include "magic_enum.h"

//stuff to be able to print type info
#include <exception>
#include <typeinfo>
#include <cxxabi.h>
#include <regex>

using namespace magic_enum::ostream_operators;
using namespace boost::program_options;


//#include "absl/container/flat_hash_map.h"

// Takes a specific numered location in the lattice and distance d and specifies the relavant stabilizer indices, type = 0 for vertices and type = 1 for plaquettes

using namespace CodeTools;

#ifdef TI_OFF
constexpr bool ti_enabled = false;
#else
constexpr bool ti_enabled = true;
#endif

constexpr int distance = DISTANCE;
int limit = 100;
std::vector<double> p_series{};
constexpr auto initial_state = single_pauli_state::ZERO;
int nc = 40;
using PlanarBase = Planar::PlanarBase<distance>;
using SimBase = SIM<PlanarBase::TOTAL_QUBITS>;
// using SimBase = ChpSimulator<PlanarBase::TOTAL_QUBITS>;
using Sim = AutoBalancedErrorSim<SimBase>;
using Planar_T = Planar::Planar<distance, Sim>;



auto intRand() {
    static thread_local std::random_device* generator = nullptr;
    if (!generator){
        generator = new std::random_device("/dev/urandom");
    }
    return (*generator)();
}

std::shared_ptr<std::mt19937_64> getRandGen(){
    std::seed_seq seq{intRand(), intRand(), intRand(), intRand()};
    std::shared_ptr<std::mt19937_64> gen = std::make_shared<std::mt19937_64>(seq);
    return gen;
}

std::string getSimType(){
    int status;
    char* realname;
    auto rand_gen = getRandGen();
    SimBase simulator_base{rand_gen};
    const std::type_info& ti = typeid(simulator_base);
    realname = abi::__cxa_demangle(ti.name(), 0, 0, &status);
    std::string realname_s = std::string(realname);
    free(realname);

    std::regex filt_re("std::mersenne_twister_engine<[^>]*>");
    std::string retval = std::regex_replace(realname_s, filt_re, "std::mt19937_64");
    return retval;
}

Planar_T get_sim(double p, single_pauli_state st){
    //s << "Setting up state" << std::endl;
    auto rand_gen = getRandGen();
    SimBase simulator_base{rand_gen};
    Sim simulator{std::move(simulator_base), rand_gen, p};
    //Planar_T surface_code(simulator, single_pauli_state::ZERO);
    return Planar_T(std::move(simulator), st, p);
}

Planar_T get_sim(double p, single_pauli_state st,  double alpha, double beta, double gamma){
    //s << "Setting up state" << std::endl;
    auto rand_gen = getRandGen();
    SimBase simulator_base{rand_gen};
    Sim simulator{std::move(simulator_base), rand_gen, p};
    //Planar_T surface_code(simulator, single_pauli_state::ZERO);
    return Planar_T(simulator, st, alpha, beta, gamma, p);
}

auto get_sim_prepared(double p, std::ostream& s, single_pauli_state st, Planar_T::ErrorSyndrome measurement, double alpha, double beta, double gamma){
    // //s << "Setting up state" << std::endl;
    // Sim simulator{getRandGen(), p, 0};
    // //Planar_T surface_code(simulator, single_pauli_state::ZERO);
    // Planar_T surface_code(simulator, single_pauli_state::ZERO, 0.434, 0.434, 0.0);
    Planar_T surface_code(get_sim(p,st,  alpha, beta, gamma));
    //std::cout << std::string(surface_code.sim) << std::endl;
    double measurement_prob = 1.0;
    Planar_T::T_2DMatch x_match{}, z_match{};

    if(ti_enabled){
        measurement_prob = surface_code.perform_round_set(measurement);
        //std::cout << std::string(surface_code.sim) << std::endl;

        x_match = surface_code.match_planar_2D(measurement, 0);
        z_match = surface_code.match_planar_2D(measurement, 1);
        // surface_code.applyMatch2D(x_match, 0);
        // surface_code.applyMatch2D(z_match, 1);
    } else {
        measurement_prob = surface_code.perform_round_set(measurement);

        x_match = surface_code.match_planar_2D(measurement, 0);
        z_match = surface_code.match_planar_2D(measurement, 1);

        surface_code.applyMatch2D(x_match, 0);
        surface_code.applyMatch2D(z_match, 1);
    }



    //surface_code.simple_correct(measurement);
    //s << measurement_prob << " " << measurement << std::endl;
    return std::make_tuple(measurement_prob, surface_code , x_match, z_match);
}

auto get_sim_prepared(double p, std::ostream& s, single_pauli_state st, uint64_t measurement_id, double alpha, double beta, double gamma){
    auto measurement = Planar_T::getSyndromeFromInt(measurement_id);
    return get_sim_prepared(p,s,st,measurement, alpha, beta, gamma);
}

Planar_T::ErrorSyndrome get_syndrome_from_match(const Planar_T::T_2DMatch & input_match_x, const Planar_T::T_2DMatch & input_match_z){
    Planar_T::ErrorSyndrome s{};

    for (auto m : input_match_x){
        for(auto pt : m){
            if(pt[0] >= 0 && pt[1] >= 0 && pt[0] < Planar_T::LATTICE_DIMENSION  && pt[1] < Planar_T::LATTICE_DIMENSION){
                s[0].push_back(pt);
            }
        }
    }

    for (auto m : input_match_z){
        for(auto pt : m){
            if(pt[0] >= 0 && pt[1] >= 0 && pt[0] < Planar_T::LATTICE_DIMENSION  && pt[1] < Planar_T::LATTICE_DIMENSION){
                s[1].push_back(pt);
            }
        }
    }
    return s;
}

Planar_T::T_2DMatch get_base_match(Planar_T::T_3DMatch input_match){
    Planar_T::T_2DMatch base_m{};

    for(auto m: input_match){
        if((m[0][2] == 0) && (m[1][2] == 0)){
            std::array<int,2> inner_first({m[0][0], m[0][1]});
            std::array<int,2> inner_second({m[1][0], m[1][1]});
            if(inner_first != inner_second){
                base_m.push_back({inner_first, inner_second});
            }
        }
    }
    return base_m;
}

Planar_T::T_3DMatch get_non_base_match(const Planar_T::T_3DMatch &input_match){
    Planar_T::T_3DMatch base_m{};

    for(auto m: input_match){
        if((m[0][2] != 0) || (m[1][2] != 0)){
            // if(m[0][2] == 0){
            //     std::array<int,3> inner_first({m[1][0], m[1][1], 0});
            //     std::array<int,3> inner_second({m[1][0], m[1][1], m[1][2]});
            //     base_m.push_back({inner_first, inner_second});
            // } else if (m[1][2] == 0){
            //     std::array<int,3> inner_first({m[0][0], m[0][1], m[0][2]});
            //     std::array<int,3> inner_second({m[0][0], m[0][1], 0});
            //     base_m.push_back({inner_first, inner_second});
            // } else{
                std::array<int,3> inner_first({m[0][0], m[0][1], m[0][2]});
                std::array<int,3> inner_second({m[1][0], m[1][1], m[1][2]});
                base_m.push_back({inner_first, inner_second});
            // }
        }
    }
    return base_m;
}

auto get_random_sim_ti(double p, std::ostream& s, single_pauli_state st, double alpha, double beta, double gamma){
    Planar_T surface_code(get_sim(p,st,alpha,beta, gamma));
    auto error_syndrome = surface_code.extractFTSyndrome();
    Planar_T::ErrorSyndrome base_syndrome{};

    Planar_T::T_2DMatch base_x{}, base_z{};

    if(ti_enabled){
        auto x_match = surface_code.match_planar_3D(error_syndrome, 0);
        auto z_match = surface_code.match_planar_3D(error_syndrome, 1);

        base_x = get_base_match(x_match);
        base_z = get_base_match(z_match);

        auto non_base_x = get_non_base_match(x_match);
        auto non_base_z = get_non_base_match(z_match);

        surface_code.applyMatch3D(non_base_x, 0);
        surface_code.applyMatch3D(non_base_z, 1);
        // std::cerr << "a" << std::flush;

        base_syndrome = get_syndrome_from_match(base_x, base_z);
    } else {
        auto x_match = surface_code.match_planar_3D(error_syndrome, 0);
        auto z_match = surface_code.match_planar_3D(error_syndrome, 1);

        surface_code.applyMatch3D(x_match, 0);
        surface_code.applyMatch3D(z_match, 1);
    }



    // s << "base syndrome: " << base_syndrome << std::endl;

    // surface_code.applyMatch2D(base_x, 0);
    // surface_code.applyMatch2D(base_z, 1);

    // base_syndrome = Planar_T::ErrorSyndrome{};

    return std::make_tuple(surface_code, base_syndrome, error_syndrome, base_x, base_z);
}


auto get_random_sim_std(double p, std::ostream& s, single_pauli_state st){
    Planar_T surface_code = get_sim(p, st);

    s << "Extracting Syndrome: ";
    auto syndrome = surface_code.extractFTSyndrome();
    s << ", " << syndrome;
    s << ", performing correction" << ", ";
    surface_code.correctFTSyndrome(syndrome);
    Planar_T::ErrorSyndrome base_syndrome{};

    return std::make_tuple(surface_code, base_syndrome);
}


std::vector<std::tuple<double, Planar_T, Planar_T::ErrorSyndrome>> get_list_of_sims(double p, std::ostream& s,single_pauli_state st, double alpha, double beta, double gamma){
    auto s_f = get_sim(0,st);
    //s << std::string(s_f.sim) << std::endl;

    std::vector<std::tuple<double, Planar_T, Planar_T::ErrorSyndrome>>  sim_list{};
    for(int64_t i = 0; i < (1<<Planar_T::REQUIRED_ANCS); ++i){
        auto [measurement_prob, sim, mx, mz] = get_sim_prepared(p,s, st,i, alpha, beta, gamma);
        auto meas = Planar_T::getSyndromeFromInt(i);
        //double measurement_prob = 1;
        //auto sim = get_sim(p,s);
        //auto meas = sim.perform_round();
        if(measurement_prob > 1E-14){
            auto s_c = sim;
            std::cout << measurement_prob << " " << meas << std::endl;
            std::cout << std::string(sim.sim) << std::endl;
            sim_list.emplace_back(measurement_prob, sim, meas);
            //std::vector<std::tuple<double, Planar_T::ErrorSyndrome, Planar_T>> sims;
        }
    }
    return sim_list;
}

void get_ti_probs(single_pauli_state st, double alpha, double beta, double gamma){
    std::vector<double> results{};
    std::vector<int> counts{};
    std::vector<double> tallies{};
    std::vector<double> times{};
    auto total_start = std::chrono::system_clock::now();
    std::cout << "    \"progress\" : [\n        ";
    for (double p : p_series) {
        auto start_time = std::chrono::system_clock::now();
        std::cout <<"\"" << p << std::flush;
        std::vector<std::pair<double, bool>> fidelities(limit, {0.0, false});
        #pragma omp parallel for schedule(dynamic, 1) firstprivate(p, limit, st, alpha, beta, gamma) num_threads(nc)
        for (int i = 0; i < limit; i++) {
            std::stringstream obuffer = std::stringstream();
            std::stringstream buffer = std::stringstream();

            auto [base_surface_code, base_measurement, orig_measurement, tst_x, tst_z] = get_random_sim_ti(p,std::cout, st, alpha, beta, gamma);
            auto [p_meas, reference_surface_code, ref_x, ref_z] = get_sim_prepared(0, buffer, st, base_measurement ,alpha, beta, gamma);

            // base_surface_code.sim.get_error_free().normalize();
            // std::cout << std::string(base_surface_code.sim) << std::endl;
            double ip = base_surface_code.sim.getInnerProduct(reference_surface_code.sim);
            buffer << "Reference Probability: " << p_meas << std::endl;
            buffer << "Inner product to reference state: " << ip << std::endl;

            // std::cout << buffer.str() << std::endl;

            // buffer << "Fidelity for test i: ";
            // double ip = surface_code_sim1.sim.getInnerProduct(surface_code_sim2.sim);
            double F = sqrt(ip);
            //std::cerr << F << std::endl;
            if (!std::isnan(F) && (p_meas > 1E-12)) {
                fidelities[i] = std::make_pair(F, true);
                if ((1-F) > 1E-10 && p == 0.0){
                    std::cerr << '\n' << 1-F << '\n' << base_measurement << '\n' <<  orig_measurement 
                              << '\n' << tst_x  <<"  "<< tst_z << '\n' << ref_x << ' ' << ref_z << std::endl;
                }
            }
            if(limit < 20 || (i % (limit/20)) == 0) {
                double tally = 0.0;
                int count = 0;
                for (const auto&[f, set]: fidelities) {
                    if (set) {
                        count++;
                        tally += f;
                    }
                }
                std:: stringstream b{};
                b << p << " Test " << i << " Running estimate: " << 1 - tally / double(count);
                // std::puts(b.str().c_str());
                std::cout << "." << std::flush;
            }
            
        }

        double tally = 0.0;
        int count = 0;
        for (const auto&[f, set]: fidelities) {
            if (set) {
                count++;
                tally += f;
            }
        }
        std::cout << ' ' << 1 - (tally / count) << ' ' << count << "\",\n        " <<  std::flush;
        results.push_back(1 - (tally / count));
        counts.push_back(count);
        tallies.push_back(tally);
        auto end_time = std::chrono::system_clock::now();
        std::chrono::duration<double> duration = end_time - start_time;
        times.push_back(duration.count());
    }
    auto total_end = std::chrono::system_clock::now();
    std::chrono::duration<double> total_time = total_end - total_start;

    std::cout << "\"\"],\n    \"config\": {\n        \"Distance\" : "<< distance << ",\n"
              << "        \"samples\" : " << limit << ",\n"
              << "        \"alpha\" : " << alpha << ",\n"
              << "        \"beta\" : " << beta  << ",\n"
              << "        \"gamma\" : " << gamma << ",\n"
              << "        \"sim_type\" : \"" << getSimType() << "\", \n"
              << "        \"ti_on\" : " << (ti_enabled ? "true": "false") << ",\n"
              << "        \"state\" : \"" << st << "\"\n    },\n";
    std::cout << "    \"physicalRates\" : "<< p_series 
           << ",\n    \"logicalRates\" : "<< results 
           << ",\n    \"tallies\" : " << tallies
           << ",\n    \"counts\" : " << counts 
           << ",\n    \"times\" : " << times
           << ",\n    \"totalTime\" : " << total_time.count()
           << "\n" << std::flush;


}


// void print_logical_probs(Planar_T sim){
//     std::cout << "Probability of a parity 1 logical operator Z";
//     std::vector<int> z_operator(distance);
//     for (int i = 0; i < distance; i++) {
//         z_operator[i] = i * 2;
//     }
//     std::cout << z_operator << ": ";
//     std::cout << sim.sim.getParityOneMeasurementProbability(z_operator) << std::endl;

//     std::cout << "Probability of a parity 1 logical operator X";
//     std::vector<int> x_operator(distance);
//     for (int i = 0; i <= distance; i++) {
//         x_operator[i] = i * PlanarBase::LATTICE_DIMENSION * 2;
//     }


//     // sim.sim.enable_errors = 0;
//     // for (int q : x_operator) {
//     //     sim.sim.h(q);
//     // }
//     // std::cout << x_operator << ": ";
//     // std::cout << sim.sim.getParityOneMeasurementProbability(x_operator) << std::endl;
//     // for (int q : x_operator) {
//     //     sim.sim.h(q);
//     // }
//     // sim.sim.enable_errors = 1;
// }

// void get_probs(single_pauli_state st){
//     std::vector<double> results{};
//     Planar_T surface_code_sim1 = get_sim(0.0, st);

//     //Extract 2D syndrome
//     Planar_T::ErrorSyndrome syndrome = surface_code_sim1.perform_round();
//     std::cout << "Measured Syndrome " << syndrome;
//     std::vector<int> syndrome_indexes{};
//     for (const auto &stab_type_list : syndrome) {
//         for (const auto &axion : stab_type_list) {
//             int index = (axion[0] * Planar_T::LATTICE_DIMENSION + axion[1]) / 2;
//             syndrome_indexes.push_back(index);
//         }
//     }

//     std::cout << "Syndrome Indexes" << syndrome_indexes << std::endl;
//     std::cout << "Performing correction" << std::endl;
//     surface_code_sim1.correctFTSyndrome(surface_code_sim1.extractFTSyndrome());


//     std::cout << "Confirm Stabilisers" << std::endl;
//     auto confirmation = surface_code_sim1.perform_round();
//     std::cout << confirmation << std::endl;

//     print_logical_probs(surface_code_sim1);

//     for (double p : p_series) {
//         std::vector<std::pair<double, bool>> fidelities(limit, {0.0, false});
//         #pragma omp parallel for schedule(dynamic, 1) firstprivate(surface_code_sim1, p, limit, st) num_threads(nc)
//         for (int i = 0; i < limit; i++) {
//             std::stringstream obuffer = std::stringstream();
//             std::stringstream buffer = std::stringstream();

//             Planar_T surface_code_sim2 = get_sim(p, st);

//             buffer << "Extracting Syndrome: ";
//             auto syndrome = surface_code_sim2.extractFTSyndrome();
//             buffer << ", " << syndrome;
//             buffer << ", performing correction" << ", ";
//             surface_code_sim2.correctFTSyndrome(syndrome);

//             //buffer << "Extracting Syndrome again: ";
//             //auto syndrome2 = surface_code_sim2.extractFTSyndrome();
//             //buffer << ", " << syndrome2;
//             //buffer << ", performing correction" << ", ";
//             //surface_code_sim2.correctFTSyndrome(syndrome2);

//             //buffer << ", Flawless 2D Correction";
//             //surface_code_sim2.sim.enable_errors = 0;
//             //auto flawless2 = surface_code_sim2.perform_round();
//             //buffer << flawless2;
//             //surface_code_sim2.simple_correct(flawless2);
//             //surface_code_sim2.sim.enable_errors = 1;

//             buffer << "Fidelity for test i: ";
//             double ip = surface_code_sim1.sim.getInnerProduct(surface_code_sim2.sim);
//             double F = sqrt(ip);
//             if (!std::isnan(F)) {
//                 fidelities[i] = std::make_pair(F, true);
//             }

//             if((i % (limit/20)) == 0) {
//                 double tally = 0.0;
//                 int count = 0;
//                 for (const auto&[f, set]: fidelities) {
//                     if (set) {
//                         count++;
//                         tally += f;
//                     }
//                 }
//                 std:: stringstream b{};
//                 b << p << " Test " << i << " Running estimate" << 1 - tally / double(count);
//                 // std::puts(b.str().c_str());
//                 std::cout << "." << std::flush;
//             }
//         }

//         double tally = 0.0;
//         int count = 0;
//         for (const auto&[f, set]: fidelities) {
//             if (set) {
//                 count++;
//                 tally += f;
//             }
//         }
//         std::cout << 1 - (tally / count) << std::endl;
//         results.push_back(1 - (tally / count));
//     }
//     std::cout << "Distance "<< distance << " code with " << limit << " samples with " << st << " state" << std::endl;
//     std::cout << "physical rates "<< p_series << "\n" << "Logical rates "<< results << std::endl;
// }

StateSimulator<1> get_transverse_state(single_pauli_state st, double alpha, double beta, double gamma)
{
    StateSimulator<1> sim(getRandGen(), 0, 0);
    switch(st){
        case single_pauli_state::ZERO:
            break;
        case single_pauli_state::ONE:
            assert(false);// not supported
            break;
        case single_pauli_state::PLUS:
            sim.h(0);
            break;
        case single_pauli_state::MINUS:
            sim.h(0);
            assert(false);
            break;
        case single_pauli_state::I_MINUS:
            assert(false);// not supported
            break;
        case single_pauli_state::I_PLUS:
            assert(false);// not supported
            break;
    }

    sim.y(0, alpha);
    sim.z(0, beta);
    sim.y(0, gamma);

    return sim;
}

/* Function used to check that 'opt1' and 'opt2' are not specified
   at the same time. */
void conflicting_options(const variables_map& vm, 
                         const char* opt1, const char* opt2)
{
    if (vm.count(opt1) && !vm[opt1].defaulted() 
        && vm.count(opt2) && !vm[opt2].defaulted())
        throw std::logic_error(std::string("Conflicting options '") 
                          + opt1 + "' and '" + opt2 + "'.");
}

/* Function used to check that of 'for_what' is specified, then
   'required_option' is specified too. */
void option_dependency(const variables_map& vm,
                        const char* for_what, const char* required_option)
{
    if (vm.count(for_what) && !vm[for_what].defaulted())
        if (vm.count(required_option) == 0 || vm[required_option].defaulted())
            throw std::logic_error(std::string("Option '") + for_what 
                              + "' requires option '" + required_option + "'.");
}

#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
int main(int argc, char* argv[]){
    // auto sim_list = get_list_of_sims(0.0,std::cout, single_pauli_state::ZERO, 0.434, 0.434, 0.0);
    // auto zero_list = get_list_of_sims(0.0,std::cout, single_pauli_state::ZERO,0.0,0.0,0.0);

    // get_random_sim(0.0,std::cout, single_pauli_state::ZERO, 0.434, 0.434, 0.0);
    // auto [base_surface_code, base_measurement] = get_random_sim(0.002,std::cout, single_pauli_state::ZERO, 0.434, 0.434, 0.0);
    // auto [p_meas, reference_surface_code] = get_sim_prepared(0, std::cout, single_pauli_state::ZERO, base_measurement,0.434, 0.434, 0.0);

    // double ip = base_surface_code.sim.getInnerProduct(reference_surface_code.sim);
    // std::cout << "Reference Probability: " << p_meas << std::endl;
    // std::cout << "Inner product to reference state: " << ip << std::endl;

    auto st = single_pauli_state::ZERO;
    std::vector<double> alphas;
    
    try {
        double alpha_start = 0.0;
        double alpha_end = 0.5;
        int a_n = 1;

        double p_start = 0.1;
        double p_end = 0.01;
        int p_n = 3;
        options_description desc("Allowed Options");
        desc.add_options()
        ("help,h", "print usage message")\
        ("alpha-start", value<double>(&alpha_start)->default_value(0.0), "start value of alpha sweep diviced by pi")
        ("alpha-end", value<double>(&alpha_end)->default_value(0.5), "end value of alpha sweep ")
        ("alpha-n", value<int>(&a_n)->default_value(1), "number of values in alpha sweep, must be greater than 1")
        ("alphas,a", value<std::vector<double> > (&alphas)->multitoken())
        ("limit,l", value<int>(&limit)->default_value(1000), "number of samples in a datapoint")
        ("probs,p", value<std::vector<double> > (&p_series)->multitoken(), "probabilities to sample")
        ("probs-start", value<double>(&p_start)->default_value(0.1), "start value of logarithmic probability sweep")
        ("probs-end", value<double>(&p_end)->default_value(0.01), "end value of logarithmic probability sweep")
        ("probs-n", value<int>(&p_n)->default_value(4), "the number of probabilities to sample in each sweep")
        ("num-cores", value<int>(&nc)->default_value(omp_get_num_procs()), "the number of cores to use")
        ;

        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);

        if(vm.count("help")){
            std::cout << desc << std::endl;
            return EXIT_SUCCESS;
        }

        conflicting_options(vm, "alphas", "alpha-start");
        conflicting_options(vm, "alphas", "alpha-end");
        conflicting_options(vm, "alphas", "alpha-n");

        notify(vm);
        if(!(vm.count("alphas"))){
            double alpha = std::numbers::pi * alpha_start;
            for (int i = 0; i < a_n ; i++){
                if(a_n > 1){
                    double sample = (alpha_end - alpha_start) * (static_cast<double>(i)) / (static_cast<double>(a_n-1)) + alpha_start;
                    alpha = std::numbers::pi * sample;
                }
                alphas.push_back(alpha);
            }
        }

        if(!(vm.count("probs"))){
            double start = std::log(p_start);
            double end = std::log(p_end);
            for (int i = 0; i < p_n ; i++){
                double p = std::exp(start);
                if(p_n > 1){
                    double sample = (end - start) * (static_cast<double>(i)) / (static_cast<double>(p_n-1)) + start;
                    p =  std::exp(sample);
                }
                p_series.push_back(p);
            }
        }
    } catch (std::exception& e){
        std::cerr << e.what() << "\n";
        return EXIT_FAILURE;
    }
    std::cout << "[\n";
    bool first = true;
    for (double alpha : alphas){
        if(first){
            first = false;
        } else{
            std::cout << ",\n";
        }
        std::cout << "  {\n";
        std::cerr << "alpha = " << alpha << std::endl;
        std::cerr <<"transverse_state: " << std::string(get_transverse_state(st, alpha, 0.0,0.0)) << std::endl;
        get_ti_probs(st, alpha, 0.0, 0.0);
        std::cout << "  }";
        // std::cout <<"transverse_state: " << std::string(get_transverse_state(st, alpha, 0.0,0.0)) << std::endl;
    }
    std::cout << "\n]\n";
}
#pragma clang diagnostic pop


/*-----------------------end-------------------------------*/

