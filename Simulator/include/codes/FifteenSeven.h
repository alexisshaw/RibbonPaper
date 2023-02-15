//
// Created by 13114347 on 10/06/2020.
//

#ifndef STEANE_SIM_FIFTEENSEVEN_H
#define STEANE_SIM_FIFTEENSEVEN_H

#include <utility>
#include <array>
#include <algorithm>
#include <map>
#include "sim/chp.h"
#include "BaseTypes.h"
#include "ErrorModels/TrivialErrorModel.h"
#include "absl/container/flat_hash_map.h"

namespace CodeTools::FifteenSeven {
    struct FifteenSevenBase {
        static constexpr int REQUIRED_QUBITS{15};
        static constexpr int REQUIRED_ANCS{3};
        static constexpr int NUM_LOGICAL{7};
        static constexpr int LEN_LOGICAL{5};
        static constexpr int NUM_STABS{4};
        static constexpr int LEN_STABS{8};

        using ErrorSyndrome = std::array<short, 3>;
        using FTErrorSyndrome = std::tuple<ErrorSyndrome, ErrorSyndrome>;
        using InitialState = std::array<single_pauli_state, NUM_LOGICAL>;

        static std::map<FTErrorSyndrome, error_correction_list> correction_map;

        static inline InitialState getProductInitialState(single_pauli_state s_basis){
            InitialState initial_state{};
            for(auto& inner_state: initial_state){
                inner_state = s_basis;
            }
            return initial_state;
        }

        // Static Members
        constexpr static std::array<int, REQUIRED_ANCS> ANC{15, 16, 17};
        constexpr static std::array<std::array<int, LEN_STABS>, NUM_STABS> STABS{{
                         {{7, 8, 9, 10, 11, 12, 13, 14}},
                         {{3, 4, 5, 6, 11, 12, 13, 14}},
                         {{1, 2, 5, 6, 9, 10, 13, 14}},
                         {{0, 2, 4, 6, 8, 10, 12, 14}}
                 }};
        constexpr static std::array<std::array<int, REQUIRED_QUBITS - LEN_STABS>, NUM_STABS> NOT_STAB{{
            {{0, 1, 2, 3, 4, 5,  6}},
            {{0, 1, 2, 7, 8, 9,  10}},
            {{0, 3, 4, 7, 8, 11, 12}},
            {{1, 3, 5, 7, 9, 11, 13}}
                                                                                                      }};
        constexpr static std::array<std::array<int, LEN_LOGICAL>, NUM_LOGICAL> QUBIT_BASIS{{
            {{0x0, 0x1, 0x3, 0x7, 0xE}},
            {{0x0, 0x1, 0x4, 0x9, 0xB}},
            {{0x0, 0x1, 0x5, 0xA, 0xD}},
            {{0x0, 0x1, 0x6, 0x8, 0xC}},
            {{0x0, 0x3, 0x5, 0x8, 0x9}},
            {{0x0, 0x3, 0x6, 0xB, 0xD}},
            {{0x0, 0x7, 0x9, 0xC, 0xD}}
                                                                                           }};
    };

    template<typename QubitSimulator=ChpSimulator<FifteenSevenBase::REQUIRED_QUBITS + FifteenSevenBase::REQUIRED_ANCS>>
    class FifteenSeven :public FifteenSevenBase{
    public:
        QubitSimulator sim{};
        const InitialState initialState;
        absl::flat_hash_map<FTErrorSyndrome, error_correction_list> my_correction_map;

    public:
        explicit FifteenSeven(const QubitSimulator& sim_in, const InitialState& state, bool do_init=true);
        explicit FifteenSeven(QubitSimulator&& sim_in, const InitialState& state, bool do_init=true);

        ErrorSyndrome perform_round();

        void simple_correct(ErrorSyndrome error_syndrome);

        FTErrorSyndrome extractFTSyndrome();
        void correctFTSyndrome(FTErrorSyndrome s);

        [[nodiscard]] bool confirm_no_error() const;

        explicit operator std::string();

        template<typename S>
        friend
        class FifteenSeven;
    private:
        measurement_state initialize();

        void initialization_rotation(InitialState s);

        void applyCorrection(const error_correction_list &correction) {
            for(auto corr_to_fix: correction){
                sim.inject_single_pauli_error(corr_to_fix.first, corr_to_fix.second);
            }
        }

        template <typename list>
        void idleExcept(list l);

        std::array<measurement_state, 3> measure_stabilizers(int stab_no);
    };

    template<typename QubitSimulator>
    inline FifteenSeven<QubitSimulator>::FifteenSeven(const QubitSimulator& sim_in, const InitialState& state, bool do_init):
    sim{sim_in}, initialState{state}, my_correction_map{correction_map.begin(), correction_map.end()}{
        auto& base = sim.get_error_free();
        if(do_init) {
            initialize();
            initialization_rotation(state);
            base.measure_and_reset_to_state(ANC[0], CodeTools::single_pauli_state::ZERO);
            base.measure_and_reset_to_state(ANC[1], CodeTools::single_pauli_state::PLUS);
            base.measure_and_reset_to_state(ANC[2], CodeTools::single_pauli_state::PLUS);
        }
    }


    template<typename QubitSimulator>
    inline FifteenSeven<QubitSimulator>::FifteenSeven(QubitSimulator &&sim_in, const InitialState& state, bool do_init):
    sim{std::move(sim_in)}, initialState{state}, my_correction_map{correction_map}{
        auto& base = sim.get_error_free();
        if(do_init) {
            initialize();
            initialization_rotation(state);
            base.measure_and_reset_to_state(ANC[0], CodeTools::single_pauli_state::ZERO);
            base.measure_and_reset_to_state(ANC[1], CodeTools::single_pauli_state::PLUS);
            base.measure_and_reset_to_state(ANC[2], CodeTools::single_pauli_state::PLUS);
        }
    }

    template<typename QubitSimulator>
    FifteenSevenBase::FTErrorSyndrome FifteenSeven<QubitSimulator>::extractFTSyndrome() {
        ErrorSyndrome syndrome1{perform_round()};
        ErrorSyndrome syndrome2{};
        if(syndrome1 != ErrorSyndrome{{0,0,0}}){
            syndrome2 = perform_round();
        }
        return std::make_tuple(syndrome1, syndrome2);
    }

    template<typename QubitSimulator>
    inline typename FifteenSeven<QubitSimulator>::ErrorSyndrome FifteenSeven<QubitSimulator>::perform_round() {
        ErrorSyndrome stab_measurements{0, 0, 0};
        for (int i = 0; i < NUM_STABS; i++) {
            std::array<measurement_state, 3> temp = measure_stabilizers(i);
            for (const int j: {0, 1, 2}) {
                stab_measurements[j] = (2u * stab_measurements[j]) + (static_cast<unsigned int>(temp[j]) & 1u);
            }
        }
        return stab_measurements;
    }

    template<typename QubitSimulator>
    inline void FifteenSeven<QubitSimulator>::simple_correct(ErrorSyndrome error_syndrome) {
        auto& base = sim.get_error_free();
        if (error_syndrome[0] != 0) {
            base.x(error_syndrome[0] - 1);
        }
        if (error_syndrome[1] != 0) {
            base.z(error_syndrome[1] - 1);
        }
    }


    template<typename QubitSimulator>
    inline std::array<measurement_state, 3>
    FifteenSeven<QubitSimulator>::measure_stabilizers(int stab_no) {
        sim.cx(STABS[stab_no][7], ANC[0]);
        sim.moment();

        sim.cx(ANC[1], STABS[stab_no][0]);
        sim.moment();

        sim.cx(ANC[2], ANC[0]);
        sim.moment();

        sim.cz(ANC[1], ANC[2]);
        sim.cx(STABS[stab_no][5], ANC[0]);
        sim.moment();

        sim.cx(ANC[1], STABS[stab_no][1]);
        sim.moment();

        sim.cx(STABS[stab_no][6], ANC[0]);
        sim.moment();

        sim.cx(ANC[1], STABS[stab_no][4]);
        sim.cx(STABS[stab_no][3], ANC[0]);
        sim.moment();

        sim.cx(ANC[1], STABS[stab_no][2]);
        sim.moment();

        sim.cx(STABS[stab_no][4], ANC[0]);
        sim.moment();

        sim.cx(ANC[1], STABS[stab_no][6]);
        sim.cx(STABS[stab_no][2], ANC[0]);
        sim.moment();

        sim.cx(ANC[1], STABS[stab_no][3]);
        sim.cx(STABS[stab_no][1], ANC[0]);
        sim.moment();

        sim.cx(ANC[1], STABS[stab_no][5]);
        sim.moment();

        sim.cx(ANC[2], ANC[0]);
        sim.moment();

        sim.cz(ANC[1], ANC[2]);
        sim.cx(STABS[stab_no][0], ANC[0]);
        sim.moment();

        sim.cx(ANC[1], STABS[stab_no][7]);
        sim.moment();

        std::array<measurement_state, 3> retval{{sim.measure_in_basis_and_reset_to_state(ANC[0], single_pauli_basis::Z, single_pauli_state::ZERO),
                                                        sim.measure_in_basis_and_reset_to_state(ANC[1], single_pauli_basis::X, single_pauli_state::PLUS),
                                                        sim.measure_in_basis_and_reset_to_state(ANC[2], single_pauli_basis::X, single_pauli_state::PLUS)}};
        sim.moment();

        return retval;
    }




    template<typename QubitSimulator>
    inline measurement_state FifteenSeven<QubitSimulator>::initialize() {
        auto& base = sim.get_error_free();
        for (const auto &stab_qubits : STABS) {
            base.measure_and_reset_to_state(stab_qubits[0], CodeTools::single_pauli_state::PLUS);
            base.cx(stab_qubits[0], ANC[0]);
            for (int j = 1; j < LEN_STABS; j++) {
                base.cx(stab_qubits[0], stab_qubits[j]);
            }
        }
        for (int i = 0; i < NUM_STABS; ++i) {
            base.cx(STABS[(NUM_STABS - 1) - i][0], ANC[0]);
        }

        base.measure_and_reset_to_state(ANC[1], CodeTools::single_pauli_state::PLUS);
        base.measure_and_reset_to_state(ANC[2], CodeTools::single_pauli_state::PLUS);
        return base.measure_and_reset_to_state(ANC[2], CodeTools::single_pauli_state::ZERO);
    }

    template<typename QubitSimulator>
    inline void FifteenSeven<QubitSimulator>::initialization_rotation(InitialState s) {
        auto& base = sim.get_error_free();
        for (int logical_qubit = 0; logical_qubit < NUM_LOGICAL; ++logical_qubit) {
            auto initial_qubit_state = s[logical_qubit];
            for (int i = 0; i < 5; i++) {
                switch (initial_qubit_state) {
                    case single_pauli_state::ZERO:
                    case single_pauli_state::ONE:
                        base.cz(FifteenSeven<TrivialErrorModel>::ANC[1], QUBIT_BASIS[logical_qubit][i]);
                        break;
                    case single_pauli_state::PLUS:
                    case single_pauli_state::MINUS:
                        base.cx(FifteenSeven<TrivialErrorModel>::ANC[1], QUBIT_BASIS[logical_qubit][i]);
                        break;
                    case single_pauli_state::I_PLUS:
                    case single_pauli_state::I_MINUS:
                        base.cy(FifteenSeven<TrivialErrorModel>::ANC[1], QUBIT_BASIS[logical_qubit][i]);
                        break;
                    default:
                        // should not occur.
                        break;
                }
            }
            measurement_state val = base.measure_in_basis_and_reset_to_state(
                    FifteenSeven<TrivialErrorModel>::ANC[1],
                    CodeTools::single_pauli_basis::X,
                    CodeTools::single_pauli_state::PLUS);

            switch (initial_qubit_state) {
                case single_pauli_state::ZERO:
                    if (val == measurement_state::DETERMINISTIC_TRUE || val == measurement_state::RANDOM_TRUE) {
                        for (int i = 0; i < 5; i++) {
                            base.x(QUBIT_BASIS[logical_qubit][i]);
                        }
                        break;
                    }
                    break;
                case single_pauli_state::PLUS:
                    if (val == measurement_state::DETERMINISTIC_TRUE || val == measurement_state::RANDOM_TRUE) {
                        for (int i = 0; i < 5; i++) {
                            base.z(QUBIT_BASIS[logical_qubit][i]);
                        }
                    }
                    break;
                case single_pauli_state::I_PLUS:
                    if (val == measurement_state::DETERMINISTIC_TRUE || val == measurement_state::RANDOM_TRUE) {
                        for (int i = 0; i < 5; i++) {
                            base.x(QUBIT_BASIS[logical_qubit][i]);
                        }
                    }
                    break;
                case single_pauli_state::ONE:
                    if (val == measurement_state::DETERMINISTIC_FALSE || val == measurement_state::RANDOM_FALSE) {
                        for (int i = 0; i < 5; i++) {
                            base.x(QUBIT_BASIS[logical_qubit][i]);
                        }
                    }
                    break;
                case single_pauli_state::MINUS:
                    if (val == measurement_state::DETERMINISTIC_FALSE || val == measurement_state::RANDOM_FALSE) {
                        for (int i = 0; i < 5; i++) {
                            base.z(QUBIT_BASIS[logical_qubit][i]);
                        }
                    }
                    break;
                case single_pauli_state::I_MINUS:
                    if (val == measurement_state::DETERMINISTIC_FALSE || val == measurement_state::RANDOM_FALSE) {
                        for (int i = 0; i < 5; i++) {
                            base.x(QUBIT_BASIS[logical_qubit][i]);
                        }
                    }
                    break;
                default:
                    break;
            }
        }
    }

    template<typename QubitSimulator>
    inline bool FifteenSeven<QubitSimulator>::confirm_no_error() const {
        using BaseSim = typename QubitSimulator::error_free_sim;
        FifteenSeven<BaseSim> clone(sim.get_error_free(), initialState, false);
        ErrorSyndrome m1 = clone.perform_round();
        clone.simple_correct(m1);
        ErrorSyndrome m2 = clone.perform_round();
        clone.simple_correct(m2);
        ErrorSyndrome m3 = clone.perform_round();
        clone.simple_correct(m3);
        clone.sim.measure_and_reset(FifteenSeven<TrivialErrorModel>::ANC[1]);
        clone.sim.h(FifteenSeven<TrivialErrorModel>::ANC[1]);

        for (int logical_qubit = 0; logical_qubit < NUM_LOGICAL; ++logical_qubit) {
            auto initial_qubit_state = initialState[logical_qubit];

            for (int i = 0; i < LEN_LOGICAL; i++) {
                switch (initial_qubit_state) {
                    case single_pauli_state::ZERO:
                    case single_pauli_state::ONE:
                        clone.sim.cz(FifteenSeven<TrivialErrorModel>::ANC[1], QUBIT_BASIS[logical_qubit][i]);
                        break;
                    case single_pauli_state::PLUS:
                    case single_pauli_state::MINUS:
                        clone.sim.cx(FifteenSeven<TrivialErrorModel>::ANC[1], QUBIT_BASIS[logical_qubit][i]);
                        break;
                    case single_pauli_state::I_PLUS:
                    case single_pauli_state::I_MINUS:
                        clone.sim.cy(FifteenSeven<TrivialErrorModel>::ANC[1], QUBIT_BASIS[logical_qubit][i]);
                        break;
                    default:
                        // should not occur.
                        break;
                }
            }

            measurement_state val = clone.sim.measure_in_basis_and_reset_to_state(
                    FifteenSeven<TrivialErrorModel>::ANC[1],
                    CodeTools::single_pauli_basis::X,
                    CodeTools::single_pauli_state::PLUS);

            switch (initial_qubit_state) {
                case single_pauli_state::ZERO:
                case single_pauli_state::PLUS:
                case single_pauli_state::I_PLUS:
                    if (val != measurement_state::DETERMINISTIC_FALSE) {
                        return false;
                    } else {
                        break;
                    }
                case single_pauli_state::ONE:
                case single_pauli_state::MINUS:
                case single_pauli_state::I_MINUS:
                    if (val != measurement_state::DETERMINISTIC_TRUE) {
                        return false;
                    } else {
                        break;
                    }
                default:
                    // should not occur
                    return false;
            }
        }
        return true;
    }

    template<typename QubitSimulator>
    inline FifteenSeven<QubitSimulator>::operator std::string() {
        return std::string(sim);
    }

    template<typename QubitSimulator>
    void FifteenSeven<QubitSimulator>::correctFTSyndrome(FifteenSevenBase::FTErrorSyndrome s) {
        const auto& const_Corrector = correction_map;

        if(const_Corrector.count(s) == 0){
            simple_correct(std::get<1>(s));
            return;
        } else {
            const auto & correction = const_Corrector.at(s);
            applyCorrection(correction);
        }

    }

    template<typename QubitSimulator>
    template<typename list>
    void FifteenSeven<QubitSimulator>::idleExcept(list l) {
        std::array<bool, REQUIRED_QUBITS + REQUIRED_ANCS> q;
        for(int i = 0; i < REQUIRED_QUBITS + REQUIRED_ANCS; ++i){
        }

    }


}

#endif //STEANE_SIM_FIFTEENSEVEN_H
