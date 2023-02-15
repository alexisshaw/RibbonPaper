//
// Created by 13114347 on 25/05/2020.
//

#ifndef STEANE_SIM_STEANE_H
#define STEANE_SIM_STEANE_H

#include <utility>
#include <array>
#include <algorithm>
#include "sim/chp.h"
#include "BaseTypes.h"

namespace CodeTools::Steane {
    class SteaneSimBase{
    public:
        static constexpr int REQUIRED_QUBITS{7};
        static constexpr int REQUIRED_ANCS{2};
        static constexpr int NUM_LOGICAL{1};
        static constexpr int NUM_STABS {3};
        static constexpr int LEN_STABS {4};

        using ErrorSyndrome = std::pair<uint8_t, uint8_t>;
        using FTErrorSyndrome = std::array<ErrorSyndrome, 3>;
        using InitialState = single_pauli_state;

        inline static InitialState getProductInitialState(single_pauli_state s_basis){
            return s_basis;
        }

        // Static Members
        static constexpr int ANC[REQUIRED_ANCS] = {7,8};
        static constexpr const int STABS[NUM_STABS][LEN_STABS] = {{3, 4, 5, 6},
                                                                  {1, 2, 5, 6},
                                                                  {0, 2, 4, 6}};
        static constexpr const int NOT_STAB[NUM_STABS][LEN_STABS] = {{0, 1, 2},
                                                                     {0, 3, 4},
                                                                     {1, 3, 5}};
    };

    template<typename QubitSimulator=ChpSimulator<SteaneSimBase::REQUIRED_QUBITS + SteaneSimBase::REQUIRED_ANCS>>
    class SteaneSimT: public SteaneSimBase {
    public:
        QubitSimulator sim;
        const InitialState initialState;


    public:
        explicit SteaneSimT(const QubitSimulator& sim_in, const InitialState& state, bool do_init=true);
        explicit SteaneSimT(QubitSimulator&& sim_in, const InitialState& state, bool do_init=true);

        ErrorSyndrome perform_round();

        void simple_correct(const ErrorSyndrome& error_syndrome);

        FTErrorSyndrome extractFTSyndrome();

        void correctFTSyndrome(FTErrorSyndrome s);

        [[nodiscard]] bool confirm_no_error() const;

        explicit operator std::string();

        template<typename S>
        friend class SteaneSimT;

    private:
        void initialize(InitialState s);
        void fault_tolerance_correct_internal(ErrorSyndrome m1, ErrorSyndrome m2);
        std::pair<measurement_state,measurement_state> measure_stabilizers(int stab_no);
    };

    //! Constructs a QubitSimulator using a given simulator, state. If do_init then we will perform initialisation on the simulator.
    template<typename QubitSimulator>
    inline SteaneSimT<QubitSimulator>::SteaneSimT(const QubitSimulator& sim_in, const InitialState& state_in, bool do_init): sim{sim_in}, initialState{state_in}{
        if(do_init) initialize(initialState);
    }

    //! Construct a QubitSimulator, moving in the state and base_sim.
    template<typename QubitSimulator>
    SteaneSimT<QubitSimulator>::SteaneSimT(QubitSimulator&& sim_in, const InitialState& state_in, bool do_init)
            : sim{std::move(sim_in)}, initialState{state_in}{
        if(do_init) initialize(initialState);
    }

    //! Get a string representation of the SteaneSimBase, for now we just print the underlying simulator
    //! TODO: return a better representation of the simulator that includes the initialState etc...
    template<typename QubitSimulator>
    inline SteaneSimT<QubitSimulator>::operator std::string() {
        return std::string(sim);
    }

    //! Initialise the quantum code into the given state.
    template<typename QubitSimulator>
    inline void SteaneSimT<QubitSimulator>::initialize(InitialState s) {
        auto& base = sim.get_error_free();
        base.get_error_free().h(1);
        base.h(2);
        base.h(3);

        base.cx(1, 0);
        base.cx(3, 5);
        base.cx(2, 6);
        base.cx(1, 4);
        base.cx(2, 0);
        base.cx(3, 6);
        base.cx(1, 5);
        base.cx(6, 4);

        base.cx(0, ANC[0]);
        base.cx(5, ANC[0]);
        base.cx(6, ANC[0]);

        base.measure_and_reset(ANC[0]);
        base.measure_and_reset_to_state(ANC[0], single_pauli_state::PLUS);

        for (int i = 0; i != 7; i++) {
            base.rotate_qubit(i, s);
        }
    }

    //! An internal function that performs the syndrome extraction of a given by the stab_no stab_no.
    template<typename QubitSimulator>
    inline std::pair<measurement_state,measurement_state> SteaneSimT<QubitSimulator>::measure_stabilizers(const int stab_no) {
        sim.cx(STABS[stab_no][2], ANC[0]);

        sim.moment();

        sim.cx(ANC[1], STABS[stab_no][0]);

        sim.moment();

        sim.cx(STABS[stab_no][0], ANC[0]);
        sim.cx(ANC[1], STABS[stab_no][2]);

        sim.moment();

        sim.cx(STABS[stab_no][1], ANC[0]);
        sim.cx(ANC[1], STABS[stab_no][3]);

        sim.moment();

        sim.cx(ANC[1], STABS[stab_no][1]);

        sim.moment();

        sim.cx(STABS[stab_no][3], ANC[0]);

        sim.moment();

        auto r1 = sim.measure_and_reset(ANC[0]);
        auto r2 = sim.measure_in_basis_and_reset_to_state(ANC[1],single_pauli_basis::X, single_pauli_state::PLUS);

        sim.moment();

        return std::make_pair(r1,r2);
    }

    //! Performs one round of error syndrome extraction, and returns the error syndrome.
    template<typename QubitSimulator>
    inline typename SteaneSimBase::ErrorSyndrome SteaneSimT<QubitSimulator>::perform_round() {
        ErrorSyndrome stab_measurements = std::make_pair(0, 0);
        for (int i = 0; i < NUM_STABS; i++) {
            std::pair<measurement_state,measurement_state> temp = measure_stabilizers(i);
            stab_measurements.first = (2u * stab_measurements.first) + (static_cast<unsigned int>(temp.first) & 1u);
            stab_measurements.second = (2u * stab_measurements.second) + (static_cast<unsigned int>(temp.second) & 1u);
        }
        return stab_measurements;
    }

    template<typename QubitSimulator>
    SteaneSimBase::FTErrorSyndrome SteaneSimT<QubitSimulator>::extractFTSyndrome() {
        auto a1 = perform_round();
        auto a2 = perform_round();
        if (a1 == a2){
            return {a1, a2, std::make_pair(0,0)};
        } else {
            auto a3 = perform_round();
            return {a1, a2, a3};
        }
    }

    //! Checks whether there is a non-trivially corectible error on the state encoded by the error correcting code.
    template<typename QubitSimulator>
    inline bool SteaneSimT<QubitSimulator>::confirm_no_error() const {
        using BaseSim = typename QubitSimulator::error_free_sim;
        SteaneSimT<BaseSim> clone = SteaneSimT<BaseSim>(sim.get_error_free(), initialState, false);
        ErrorSyndrome m = clone.perform_round();
        clone.simple_correct(m);

        clone.sim.h(clone.ANC[0]);
        // We care here what the exact type of an InitialState is, so we breach the veil.
        for (int i = 0; i < 7; i++) {
            switch (initialState) {
                case single_pauli_state::ZERO:
                case single_pauli_state::ONE:
                    clone.sim.cz(SteaneSimBase::ANC[0], i);
                    break;
                case single_pauli_state::PLUS:
                case single_pauli_state::MINUS:
                    clone.sim.cx(SteaneSimBase::ANC[0], i);
                    break;
                case single_pauli_state::I_PLUS:
                case single_pauli_state::I_MINUS:
                    clone.sim.cy(SteaneSimBase::ANC[0], i);
                    break;
                default:
                    // should not occur.
                    break;
            }
        }

        measurement_state val = clone.sim.measure_in_basis_and_reset_to_state(
                SteaneSimBase::ANC[0],
                CodeTools::single_pauli_basis::X,
                CodeTools::single_pauli_state::ZERO);

        switch (initialState) {
            case single_pauli_state::ZERO:
            case single_pauli_state::PLUS:
            case single_pauli_state::I_PLUS:
                return val == measurement_state::DETERMINISTIC_FALSE;
            case single_pauli_state::ONE:
            case single_pauli_state::MINUS:
            case single_pauli_state::I_MINUS:
                return val == measurement_state::DETERMINISTIC_TRUE;
            default:
                return false;
        }
    }

    //! Performs a correction for an error within the weight of the code, for a flawlessly extracted error syndrome.
    template<typename QubitSimulator>
    inline void SteaneSimT<QubitSimulator>::simple_correct(const ErrorSyndrome& error_syndrome) {
        auto& base = sim.get_error_free();
        if (error_syndrome.first != 0) {
            base.x(error_syndrome.first - 1);
        }
        if (error_syndrome.second != 0) {
            base.z(error_syndrome.second - 1);
        }
    }

    //! Performs correction of a state given that a fault has occured during the extraction of error_syndromd m1,
    //! where m2 was calculated immediately after.
    template<typename QubitSimulator>
    inline void SteaneSimT<QubitSimulator>::fault_tolerance_correct_internal(ErrorSyndrome m1, ErrorSyndrome m2) {
        auto& base = sim.get_error_free();
        if (constexpr std::array<ErrorSyndrome, 6> l_a{{{6, 0}, {6, 2}, {6, 3},
                                                     {6, 4},{6, 6},{6, 7}}};
        std::find(l_a.begin(), l_a.end(), m1) != l_a.end()) {
            if (const ErrorSyndrome r1{2, 0}; m2 == r1) {
                base.x(3);
                base.x(5);
            } else if (const ErrorSyndrome r2{2, 6}; m2 == r2) {
                base.x(3);
                base.y(5);
            } else if (const ErrorSyndrome r3{2, 7}; m2 == r3) {
                base.x(3);
                base.x(5);
                base.z(6);
            } else {
                simple_correct(m2);
            }
        } else if (constexpr std::array<ErrorSyndrome, 4> l_b{{{0, 6}, {1, 6},
                                                            {4, 6},{5, 6}}};
                std::find(l_b.begin(), l_b.end(), m1) != l_b.end()) {
            if (const ErrorSyndrome r1{0, 2}; m2 == r1) {
                base.z(3);
                base.z(5);
            } else if (const ErrorSyndrome r2{4, 2}; m2 == r2) {
                base.y(3);
                base.z(5);
            } else if (const ErrorSyndrome r3{5, 2}; m2 == r3) {
                base.z(3);
                base.z(5);
                base.x(4);
            } else {
                simple_correct(m2);
            }
        } else if (
            constexpr struct {
                const std::array<ErrorSyndrome, 2> l1 {{{1, 0}, {1, 1}}};
                const std::array<ErrorSyndrome, 3> l2 {{{4, 0}, {4, 5}, {4, 7}}};
            } l_c{};
                (std::find(l_c.l1.begin(), l_c.l1.end(), m1) != l_c.l1.end()) &&
                (std::find(l_c.l2.begin(), l_c.l2.end(), m2) != l_c.l2.end())) {
            if (const ErrorSyndrome r1{4, 0}; m2 == r1) {
                base.x(2);
                base.x(6);
            } else if (const ErrorSyndrome r2{4, 7}; m2 == r2) {
                base.x(2);
                base.y(6);
            } else if (const ErrorSyndrome r3{4, 5}; m2 == r3) {
                base.x(2);
                base.x(6);
                base.z(4);
            } else {
                simple_correct(m2);
            }
        } else if (
            constexpr struct {
                const std::array<ErrorSyndrome, 4> l1 {{{0, 2}, {1, 2},
                                                         {2,2},{3,2}}};
                const std::array<ErrorSyndrome, 3> l2 {{{0, 4}, {2, 4}, {3, 4}}} ;
            } l_d{};
                (std::find(l_d.l1.begin(), l_d.l1.end(), m1) != l_d.l1.end()) &&
                (std::find(l_d.l2.begin(), l_d.l2.end(), m2) != l_d.l2.end())) {
            if (const ErrorSyndrome r1{0, 4}; m2 == r1) {
                base.z(2);
                base.z(6);
            } else if (const ErrorSyndrome r2{3, 4}; m2 == r2) {
                base.y(2);
                base.z(6);
            } else if (const ErrorSyndrome r3{2, 4}; m2 == r3) {
                base.z(2);
                base.z(6);
                base.x(1);
            } else {
                simple_correct(m2);
            }
        } else if (
            constexpr struct{
                const std::array<ErrorSyndrome, 4> l1 {{{2, 0}, {2, 1},
                                                         {2,2},{2,3}}};
                const std::array<ErrorSyndrome, 3> l2 {{{4, 0}, {4, 6}, {4, 7}}} ;
            } l_e{};
                (std::find(l_e.l1.begin(), l_e.l1.end(), m1) != l_e.l1.end()) &&
                (std::find(l_e.l2.begin(), l_e.l2.end(), m2) != l_e.l2.end())) {
            if (const ErrorSyndrome r1{4, 0}; m2 == r1) {
                base.x(1);
                base.x(5);
            } else if (const ErrorSyndrome r2{4, 6}; m2 == r2) {
                base.x(1);
                base.y(5);
            } else if (const ErrorSyndrome r3{4, 7}; m2 == r3) {
                base.x(1);
                base.x(5);
                base.z(6);
            } else {
                simple_correct(m2);
            }
        }else if (
                constexpr struct {
                const std::array<ErrorSyndrome, 2> l1 {{{0, 1}, {1, 1}}};
                const std::array<ErrorSyndrome, 3> l2 {{{0, 4}, {1, 4}, {3, 4}}} ;
            } l_f{};
                (std::find(l_f.l1.begin(), l_f.l1.end(), m1) != l_f.l1.end()) &&
                (std::find(l_f.l2.begin(), l_f.l2.end(), m2) != l_f.l2.end())) {
            if (const ErrorSyndrome r1{0, 4}; m2 == r1) {
                base.z(2);
                base.z(6);
            } else if (const ErrorSyndrome r2{3, 4}; m2 == r2) {
                base.y(2);
                base.z(6);
            } else if (const ErrorSyndrome r3{1, 4}; m2 == r3) {
                base.z(2);
                base.z(6);
                base.x(0);
            } else {
                simple_correct(m2);
            }
        } else {
            simple_correct(m2);
        }
    }

    template<typename QubitSimulator>
    inline void SteaneSimT<QubitSimulator>::correctFTSyndrome(FTErrorSyndrome s) {
        if(s[0] == s[1]){
            simple_correct(s[0]);
        } else {
            if (s[1] == s[2]){
                fault_tolerance_correct_internal(s[0], s[1]);
            } else {
                fault_tolerance_correct_internal(s[1], s[2]);
            }
        }

    }


}
#endif //STEANE_SIM_STEANE_H
