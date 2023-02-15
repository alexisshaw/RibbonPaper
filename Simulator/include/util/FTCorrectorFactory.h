//
// Created by 13114347 on 8/10/2020.
//

#ifndef STEANE_SIM_FTCORRECTORFACTORY_H
#define STEANE_SIM_FTCORRECTORFACTORY_H

#include <iostream>
#include <chrono>
#include <vector>
#include <numeric>
#include <iterator>
#include <set>
#include <map>
#include <variant>
#include "sim/chp.h"
#include "codes/FifteenSeven.h"
#include "sim/balanced_error_sim.h"
#include "sim/tee_sim.h"
#include "sim/ErrorTracker.h"
#include "stats.hpp"
#include "magic_enum.h"
#include "sim/fixed_error_sim.h"
#include "prettyprint.h"
#include "utilityFunctions.h"

namespace CorrectorFactory {
    using CodeTools::ChpSimulator;
    using CodeTools::FixedErrorSim;
    using CodeTools::ErrorTracker;
    using CodeTools::TeeSim;
    using CodeTools::BalancedErrorSim;
    using CodeTools::single_pauli_state;

    template<typename CodeTypes, template<typename sim> typename Code>
    class FTCorrectorFactory {
        static constexpr int REQUIRED_ANCS = CodeTypes::REQUIRED_ANCS;
        static constexpr int REQUIRED_QUBITS = CodeTypes::REQUIRED_QUBITS;
        static constexpr int NUM_LOGICAL = CodeTypes::NUM_LOGICAL;
        static constexpr int LEN_LOGICAL = CodeTypes::LEN_LOGICAL;
        static constexpr int NUM_STABS = CodeTypes::NUM_STABS;
        static constexpr int LEN_STABS = CodeTypes::LEN_STABS;

        using CodeBaseT = Code<FixedErrorSim<ChpSimulator<
                REQUIRED_ANCS + REQUIRED_QUBITS>>>;
        using CodeBaseTrackingCHPSimT = Code<FixedErrorSim<TeeSim<
                ChpSimulator<REQUIRED_ANCS + REQUIRED_QUBITS>, ErrorTracker<REQUIRED_ANCS + REQUIRED_QUBITS>>>>;

        using InitialState = typename CodeTypes::InitialState;
        using ErrorSyndrome = typename CodeTypes::ErrorSyndrome;
        using FTErrorSyndrome = typename CodeTypes::FTErrorSyndrome;

//        using two_error_key = std::pair<two_qubit_pauli, int>;
//        using single_error_key = std::pair<single_qubit_pauli, int>;
//        using error_key = std::variant<single_error_key, two_error_key>;

    private:
        std::shared_ptr<std::mt19937_64> rand_gen;

    public:
        FTCorrectorFactory() {
            std::random_device rd{"/dev/random"};
            std::seed_seq seq{rd(), rd(), rd(), rd()};
            rand_gen = std::make_shared<std::mt19937_64>(seq);
        }

        std::map<FTErrorSyndrome, error_correction_list> getAndTestCorrectionMap(){
            std::cout << "Getting the number of ops in a single round of syndrome extraction" << std::endl;
            const auto [num_single, num_dual] = getStabOpCount();
            std::cout << "    single_qubit_ops: " << num_single << " two_qubit_ops: " << num_dual << '\n'<< std::endl;

            std::cout << "Extracting Error Syndromes" << std::endl;
            std::multimap<FTErrorSyndrome, error_correction_list> syndromeCorrections{};

            std::cout << "  Data Errors: " << std::flush;
            syndromeCorrections.merge(extract_data_errors());
            std::cout << "done" << std::endl;

            std::cout << "  Single-Qubit Gate Errors: " << std::flush;
            syndromeCorrections.merge(extract_single_errors());
            std::cout << "done" << std::endl;

            std::cout << "  Two-Qubit Gate Errors: " << std::flush;
            syndromeCorrections.merge(extract_dual_errors());
            std::cout << "done" << std::endl;

            std::cout << "Constructing Correction map: " << std::flush;
            auto [correctionMap, correctionObstructions] = createCorrectionMap(syndromeCorrections);
            bool correctionMapGood = corrMapGood(correctionMap, correctionObstructions);

            std::cout << (correctionMapGood? "Done" : "Error") << std::endl;

            std::cout << "Testing Correction Map" << std::endl;
            std::cout << "  Data Errors: " << std::flush;
            std::cout << (test_correct_initial_state(correctionMap)? "Done":"Error") << std::endl;

            std::cout << "  Single-Qubit Gate Errors: " << std::flush;
            std::cout << (test_single_errors(correctionMap)? "Done":"Error") << std::endl;

            std::cout << "  Two-Qubit Gate Errors: " << std::flush;
            std::cout << (test_dual_errors(correctionMap)? "Done":"Error") << std::endl;
            return correctionMap;
        }
    private:
        auto getSim(long double p){
            ChpSimulator<REQUIRED_QUBITS + REQUIRED_ANCS> chp = ChpSimulator<REQUIRED_QUBITS + REQUIRED_ANCS>(rand_gen);
            BalancedErrorSim<decltype(chp)> balSim = BalancedErrorSim(p,p,0,0,chp, rand_gen);
            return balSim;
        }

        bool corrMapGood(std::map<FTErrorSyndrome, error_correction_list> &correctionMap,
                         const std::multimap<FTErrorSyndrome, error_correction_list> &correctionObstructions);

        std::pair<std::map<FTErrorSyndrome, error_correction_list>, std::multimap<FTErrorSyndrome, error_correction_list>>
        createCorrectionMap(const std::multimap<FTErrorSyndrome, error_correction_list> &combinedCorrections);

        template <typename SimT>
        error_correction_list getCorrection(const SimT& sim);

        template <typename SimT>
        void applyCorrection(SimT& sim, const error_correction_list& correction);

        std::multimap<FTErrorSyndrome, error_correction_list> extract_data_errors();

        std::multimap<std::tuple<ErrorSyndrome, ErrorSyndrome>, error_correction_list> extract_single_errors();

        std::multimap<std::tuple<ErrorSyndrome, ErrorSyndrome>, error_correction_list> extract_dual_errors();

        template <typename SimT, typename CorrectorT>
        bool ExtractSyndromesAndCheckForErrors(SimT sim, CorrectorT corrector);;

        bool test_correct_initial_state(const std::map<FTErrorSyndrome, error_correction_list>& corrector);

        bool test_single_errors(const std::map<FTErrorSyndrome, error_correction_list>& corrector);

        bool test_dual_errors(const std::map<FTErrorSyndrome, error_correction_list>& corrector);

        std::pair<int,int> getStabOpCount();

        std::vector<CodeTools::two_qubit_pauli> getTwoErrors();

        template<typename PolicyT>
        std::map<FTErrorSyndrome, error_correction_list> consolidateCorrections(std::multimap<FTErrorSyndrome, error_correction_list> corrections, PolicyT policy);

        template<typename ErrorPauliT>
        CodeBaseTrackingCHPSimT getCodeBaseTrackingCHPSim(int gate_no, InitialState initial_state, ErrorPauliT p){
            ChpSimulator<REQUIRED_ANCS + REQUIRED_QUBITS> chp_sim(rand_gen);
            ErrorTracker<REQUIRED_ANCS + REQUIRED_QUBITS> error_tracker{};
            TeeSim<decltype(chp_sim), decltype(error_tracker)> tee_sim(chp_sim, error_tracker);

            FixedErrorSim<decltype(tee_sim)> fixed_sim(tee_sim, p, gate_no);

            Code<decltype(fixed_sim)> code_sim(fixed_sim, initial_state);
            return code_sim;
        }

        template<typename ErrorPauliT>
        CodeBaseT fifteenFixedErrSimFactory(int gate_no, InitialState initial_state, ErrorPauliT p){
            ChpSimulator<REQUIRED_ANCS + REQUIRED_QUBITS> chp_sim(rand_gen);

            FixedErrorSim<decltype(chp_sim)> fixed_sim(chp_sim, p, gate_no);

            Code<decltype(fixed_sim)> code_sim(fixed_sim, initial_state);
            return code_sim;
        }
    };

    template<typename CodeTypes, template<typename sim> typename Code>
    std::vector<CodeTools::two_qubit_pauli> FTCorrectorFactory<CodeTypes, Code>::getTwoErrors() {
        std::vector<CodeTools::two_qubit_pauli> two_errors;
        for(single_qubit_pauli p1 : {single_qubit_pauli::I, single_qubit_pauli::X, single_qubit_pauli::Z, single_qubit_pauli::Y}){
            for(single_qubit_pauli p2 : {single_qubit_pauli::I, single_qubit_pauli::X, single_qubit_pauli::Z, single_qubit_pauli::Y}){
                if(!(p1 == single_qubit_pauli::I && p2 == single_qubit_pauli::I)) {
                    two_errors.emplace_back(p1, p2); 
                }
            }
        }
        return two_errors;
    }

//    template<typename CodeTypes, template<typename sim> typename Code>
//    typename FTCorrectorFactory<CodeTypes, Code>::InitialState FTCorrectorFactory<CodeTypes, Code>::state_from_int(int_fast32_t in) {
//        InitialState outState{};
//        for (int i = 0; i < 7; i++){
//            outState[i] = single_pauli_state(in % 6);
//            in = in /6;
//        }
//        return outState;
//    }

    template<typename CodeTypes, template<typename sim> typename Code>
    std::pair<int, int> FTCorrectorFactory<CodeTypes, Code>::getStabOpCount() {
        int num_single;
        int num_dual;
        ChpSimulator<REQUIRED_ANCS + REQUIRED_QUBITS> chp_sim(rand_gen);
        FixedErrorSim<decltype(chp_sim)> fixed_sim(chp_sim);
        InitialState state = CodeTypes::getProductInitialState(single_pauli_state::ZERO);
        Code<decltype(fixed_sim)> fifteen_sim(fixed_sim, state);
        fifteen_sim.perform_round();
        num_dual = fifteen_sim.sim.getDoubleCount();
        num_single = fifteen_sim.sim.getSingleCount();

        return std::make_pair(num_single, num_dual);
    }

    template<typename CodeTypes, template<typename sim> typename Code>
    bool FTCorrectorFactory<CodeTypes, Code>::test_dual_errors(
            const std::map<FTErrorSyndrome, error_correction_list> &corrector) {
        const auto [num_single, num_dual] = getStabOpCount();

        for(auto p : getTwoErrors()) {
            for (single_pauli_state s_basis: {single_pauli_state::ZERO, single_pauli_state::PLUS, single_pauli_state::I_PLUS}) {
                for(int i = 0; i < num_dual; i++){
                    //Get the Simulator
                    auto fifteen_sim = fifteenFixedErrSimFactory(i, CodeTypes::getProductInitialState(s_basis), p);

                    if(ExtractSyndromesAndCheckForErrors(fifteen_sim, corrector))
                        return false;
                }
            }
        }

        return true;
    }

    template<typename CodeTypes, template<typename sim> typename Code>
    bool FTCorrectorFactory<CodeTypes, Code>::test_single_errors(
            const std::map<FTErrorSyndrome, error_correction_list> &corrector) {
        const auto [num_single, num_dual] = getStabOpCount();

        for(single_qubit_pauli p : {single_qubit_pauli::X, single_qubit_pauli::Z, single_qubit_pauli::Y}){
            for(single_pauli_state s_basis: {single_pauli_state::ZERO, single_pauli_state::PLUS, single_pauli_state::I_PLUS}) {
                for(int i = 0; i < num_single; ++i){
                    //Get the Simulator
                    auto fifteen_sim = fifteenFixedErrSimFactory(i, CodeTypes::getProductInitialState(s_basis), p);

                    if(ExtractSyndromesAndCheckForErrors(fifteen_sim, corrector))
                        return false;
                }
            }
        }
        return true;
    }

    template<typename CodeTypes, template<typename sim>  typename Code>
    bool FTCorrectorFactory<CodeTypes, Code>::test_correct_initial_state(
            const std::map<FTErrorSyndrome, error_correction_list> &corrector) {
        for(single_qubit_pauli p : {single_qubit_pauli::X, single_qubit_pauli::Z, single_qubit_pauli::Y}) {
            for(single_pauli_state s_basis: {single_pauli_state::ZERO, single_pauli_state::PLUS, single_pauli_state::I_PLUS}) {
                for (int i = 0; i < (REQUIRED_QUBITS + REQUIRED_ANCS); ++i) {
                    //Get a new simulator - does not inject an error
                    auto fifteen_sim = fifteenFixedErrSimFactory(0, CodeTypes::getProductInitialState(s_basis),
                                                                 single_qubit_pauli::I);
                    //Manually inject an error into the qubit
                    fifteen_sim.sim.inject_single_pauli_error(i, p);
                    if(ExtractSyndromesAndCheckForErrors(fifteen_sim, corrector))
                        return false;
                }
            }
        }
        return true;
    }

    template<typename CodeTypes, template<typename sim>  typename  Code>
    template<typename SimT, typename CorrectorT>
    bool FTCorrectorFactory<CodeTypes, Code>::ExtractSyndromesAndCheckForErrors(SimT sim, CorrectorT corrector) {
        // Extract the syndrome
        FTErrorSyndrome syndrome(sim.extractFTSyndrome());
        const auto& correction = corrector.at(syndrome);

        // Apply the correction
        applyCorrection(sim, correction);

        // Apply another round of syndrome extraction and correction, this time without errors,
        // For those weight 1 ancilla initialisation errors, that propagate one round.

        // Extract the syndrome
        FTErrorSyndrome syndrome2(sim.extractFTSyndrome());
        const auto& correction2 = corrector.at(syndrome2);

        // Apply the correction
        applyCorrection(sim, correction2);

        return !sim.confirm_no_error();
    }

    template<typename CodeTypes, template<typename sim>  typename Code>
    std::multimap<std::tuple<typename FTCorrectorFactory<CodeTypes, Code>::ErrorSyndrome,
                             typename FTCorrectorFactory<CodeTypes, Code>::ErrorSyndrome>,
                             error_correction_list>
    FTCorrectorFactory<CodeTypes, Code>::extract_dual_errors() {
        const auto [num_single, num_dual] = getStabOpCount();
        std::multimap<std::tuple<ErrorSyndrome, ErrorSyndrome>, error_correction_list> corrections;

        for(auto p : getTwoErrors()) {
            for (single_pauli_state s_basis: {single_pauli_state::ZERO, single_pauli_state::PLUS, single_pauli_state::I_PLUS}) {
                for(int i = 0; i < num_dual; i++){
                    //Get the Simulator
                    auto fifteen_sim = getCodeBaseTrackingCHPSim(i, CodeTypes::getProductInitialState(s_basis), p);

                    FTErrorSyndrome syndrome(fifteen_sim.extractFTSyndrome());

                    // Insert the correction into the corrections multimap
                    corrections.insert(std::make_pair(syndrome, getCorrection(fifteen_sim)));
                }
            }
        }

        return corrections;
    }

    template<typename CodeTypes,template<typename sim> typename Code>
    std::multimap<std::tuple<typename FTCorrectorFactory<CodeTypes, Code>::ErrorSyndrome, typename FTCorrectorFactory<CodeTypes, Code>::ErrorSyndrome>, error_correction_list>
    FTCorrectorFactory<CodeTypes, Code>::extract_single_errors() {
        const auto [num_single, num_dual] = getStabOpCount();
        std::multimap<std::tuple<ErrorSyndrome, ErrorSyndrome>, error_correction_list> corrections;

        for(single_qubit_pauli p : {single_qubit_pauli::X, single_qubit_pauli::Z, single_qubit_pauli::Y}){
            for(single_pauli_state s_basis: {single_pauli_state::ZERO, single_pauli_state::PLUS, single_pauli_state::I_PLUS}) {
                for(int i = 0; i < num_single; ++i){
                    //Get the Simulator
                    auto fifteen_sim = getCodeBaseTrackingCHPSim(i, CodeTypes::getProductInitialState(s_basis), p);

                    FTErrorSyndrome syndrome(fifteen_sim.extractFTSyndrome());

                    // Insert the correction into the corrections multimap
                    corrections.insert(std::make_pair(syndrome, getCorrection(fifteen_sim)));
                }
            }
        }
        return corrections;
    }

    template<typename CodeTypes,template<typename sim> typename Code>
    std::multimap<typename FTCorrectorFactory<CodeTypes, Code>::FTErrorSyndrome, error_correction_list> FTCorrectorFactory<CodeTypes, Code>::extract_data_errors() {
        std::multimap<FTErrorSyndrome, error_correction_list> corrections;

        for(single_qubit_pauli p : {single_qubit_pauli::X, single_qubit_pauli::Z, single_qubit_pauli::Y}) {
            for(single_pauli_state s_basis: {single_pauli_state::ZERO, single_pauli_state::PLUS, single_pauli_state::I_PLUS}) {
                for (int i = 0; i < (REQUIRED_QUBITS + REQUIRED_ANCS); ++i) {
                    //Get a new simulator - does not inject an error
                    auto fifteen_sim = getCodeBaseTrackingCHPSim(0, CodeTypes::getProductInitialState(s_basis),
                                                                 single_qubit_pauli::I);
                    //Manually inject an error into the qubit
                    fifteen_sim.sim.inject_single_pauli_error(i, p);

                    // Extract the syndrome
                    FTErrorSyndrome syndrome(fifteen_sim.extractFTSyndrome());

                    // Insert the correction into the corrections multimap
                    corrections.insert(std::make_pair(syndrome, getCorrection(fifteen_sim)));
                }
            }
        }
        return corrections;
    }

    template<typename CodeTypes, template<typename sim> typename Code>
    template<typename SimT>
    void FTCorrectorFactory<CodeTypes, Code>::applyCorrection(SimT &sim, const error_correction_list &correction) {
        for(auto corr_to_fix: correction){
            sim.sim.inject_single_pauli_error(corr_to_fix.first, corr_to_fix.second);
        }
    }

    template<typename CodeTypes, template<typename sim> typename Code>
    template<typename SimT>
    error_correction_list FTCorrectorFactory<CodeTypes, Code>::getCorrection(const SimT &sim) {
        const auto& tracker = sim.sim.baseSim.altSim;
        return tracker.getCorrection();
    }



    template<typename CodeTypes, template<typename sim> typename Code>
    std::pair<std::map<typename FTCorrectorFactory<CodeTypes, Code>::FTErrorSyndrome, error_correction_list>,
              std::multimap<typename FTCorrectorFactory<CodeTypes, Code>::FTErrorSyndrome, error_correction_list>>
    FTCorrectorFactory<CodeTypes, Code>::createCorrectionMap(
            const std::multimap<FTErrorSyndrome, error_correction_list> &combinedCorrections) {
        std::multimap<FTErrorSyndrome, error_correction_list> badCorrectionAlternatives{};

        auto collisionPolicy = [&badCorrectionAlternatives](auto key, auto value, auto& map_value){
            if(std::get<0>(key) == ErrorSyndrome{{0,0,0}}){
                if(map_value.size() > value.size()){
                    badCorrectionAlternatives.insert(std::make_pair(key, map_value));
                    map_value = value;
                } else {
                    badCorrectionAlternatives.insert(std::make_pair(key, value));
                }
            } else {
                badCorrectionAlternatives.insert(std::make_pair(key, value));
            }
        };

        auto combinedConsolidatedCorrections(consolidateCorrections(combinedCorrections, collisionPolicy));
        for(auto& [syndrome, correction]: combinedConsolidatedCorrections){
            correction = simplifyErrorSyndrome<CodeBaseT>(correction);
        }


        size_t max_correction_length = 0;
        badCorrectionAlternatives = simplify_multimap(badCorrectionAlternatives);
        for(auto& [syndrome, correction]: badCorrectionAlternatives){
            correction = simplifyErrorSyndrome<CodeBaseT>(correction);
            if(correction.size() > max_correction_length)
                max_correction_length = correction.size();
        }
        badCorrectionAlternatives = simplify_multimap(badCorrectionAlternatives);

        return std::make_pair(combinedConsolidatedCorrections, badCorrectionAlternatives);
    }

    template<typename CodeTypes, template<typename sim> typename Code>
    bool FTCorrectorFactory<CodeTypes, Code>::corrMapGood(std::map<FTErrorSyndrome, error_correction_list> &correctionMap,
                                                          const std::multimap<FTErrorSyndrome, error_correction_list> &correctionObstructions) {
        bool correctionMapGood = true;
        correctionMapGood = correctionMap.count(FTErrorSyndrome()) != 0 && correctionMap.at(FTErrorSyndrome()).empty();
        for(const auto& [key, value] : correctionObstructions){
            correctionMapGood &= (key == FTErrorSyndrome() ) && value.size() == 1;
        }
        return correctionMapGood;
    }

    template<typename CodeTypes,template<typename sim> typename Code>
    template<typename PolicyT>
    std::map<typename FTCorrectorFactory<CodeTypes, Code>::FTErrorSyndrome, error_correction_list> FTCorrectorFactory<CodeTypes, Code>::consolidateCorrections(
            std::multimap<typename FTCorrectorFactory<CodeTypes, Code>::FTErrorSyndrome, error_correction_list> corrections, PolicyT policy) {
        std::map<FTErrorSyndrome, error_correction_list> consolidatedCorrections;

        for(const auto& [key, value] : corrections){
            //std::cout << key << "\t" << value.second << "\t" << value.first << std::endl;
            if(consolidatedCorrections.count(key) == 0){
                consolidatedCorrections.insert(std::make_pair(key, value));
            } else {
                if(compare_ErrorSyndrome<CodeBaseT>(consolidatedCorrections[key], value)){
                    if(value.size() < consolidatedCorrections[key].size()){
                        consolidatedCorrections.at(key) = value;
                    }
                    continue;
                } else {
                    policy(key, value, consolidatedCorrections[key]);

                }
            }
        }
        return consolidatedCorrections;
    }

};

#endif //STEANE_SIM_FTCORRECTORFACTORY_H
