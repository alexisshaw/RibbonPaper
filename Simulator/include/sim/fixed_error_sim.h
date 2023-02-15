//
// Created by 13114347 on 29/06/2020.
//

#ifndef STEANE_SIM_FIXED_ERROR_SIM_H
#define STEANE_SIM_FIXED_ERROR_SIM_H

#include <string>
#include <random>
#include "BaseTypes.h"
#include "sim/balanced_error_sim.h"

namespace CodeTools {
    class FixedErrorSimBase{

    };

    template<typename BaseSim>
    class FixedErrorSim {
    public:
        BaseSim baseSim;
    private:
        single_qubit_pauli single_error = single_qubit_pauli::I;
        two_qubit_pauli dual_error = std::make_pair(single_qubit_pauli::I, single_qubit_pauli::I);
        int single_count = 0;
        int dual_count = 0;
        int single_count_target = -1;
        int dual_count_target = -1;

        void maybe_single_error(int qubit) {
            if (single_count_target == (single_count++)) {
                baseSim.inject_single_pauli_error(qubit, single_error);
            }
        }

        void maybe_dual_error(int q1, int q2) {
            if (dual_count_target == (dual_count++)) {
                baseSim.inject_single_pauli_error(q1, dual_error.first);
                baseSim.inject_single_pauli_error(q2, dual_error.second);
            }
        }

        void maybe_measure_error_pre(int qubit) {
            if (single_count_target == (single_count++)) {
                baseSim.inject_single_pauli_error(qubit, single_error);
            }
        }

        void maybe_measure_error_post(int qubit) {
            if (single_count_target == (single_count++)) {
                baseSim.inject_single_pauli_error(qubit, single_error);
            }
        }
    public:
        static constexpr int N{BaseSim::N};

        //! A base simulator class, operations on which occur error free
        using error_free_sim = typename BaseSim::error_free_sim;

        //! This is the Default constructor, uses 4 calls to std::random to initialize a seed sequence with std:seed_seq.
        explicit FixedErrorSim(const BaseSim &b) : baseSim{b} {};

        FixedErrorSim(const BaseSim &b, single_qubit_pauli p, int target) :
                baseSim{b}, single_error(p), single_count_target(target) {};

        FixedErrorSim(const BaseSim &b, two_qubit_pauli p, int target) :
                baseSim{b}, dual_error(p), dual_count_target(target) {};
        int getSingleCount(){
            return single_count;
        }
        int getDoubleCount(){
            return dual_count;
        }

        void inject_single_pauli_error(int qubit, single_qubit_pauli p) {
            baseSim.inject_single_pauli_error(qubit, p);
        };

        //! Enables an error model sim to begin a moment;
        //! Does nothing here
        void moment(){};

        //! Performs a controlled-X gate, (CNOT gate).
        //! \param control The control qubit id.
        //! \param target The target qubit id.
        void cx(int control, int target) {
            maybe_dual_error(control, target);
            baseSim.cx(control, target);
        }

        //! Performs a controlled-Y gate
        //! \param control The control qubit id.
        //! \param target The target qubit id.
        void cy(int control, int target) {
            maybe_dual_error(control, target);
            baseSim.cy(control, target);
        }

        //! Performs a controlled-Z gate
        //! \param control The control qubit id.
        //! \param target The target qubit id.
        void cz(int control, int target) {
            maybe_dual_error(control, target);
            baseSim.cz(control, target);
        }

        //! Performs a Hadamard Gate
        //! \param qubit the qubit to perform the gate on.
        void h(int qubit) {
            maybe_single_error(qubit);
            baseSim.h(qubit);
        }

        //! Performs a Phase Gate
        //! \param qubit the qubit to perform the gate on.
        void s(int qubit) {
            maybe_single_error(qubit);
            baseSim.s(qubit);
        }

        //! Performs an inverse Phase Gate
        //! \param qubit the qubit to perform the gate on.
        void sdg(int qubit) {
            maybe_single_error(qubit);
            baseSim.sdg(qubit);
        }

        //! Performs an X-Pauli gate
        //! \param qubit the qubit to perform the gate on.
        void x(int qubit) {
            maybe_single_error(qubit);
            baseSim.x(qubit);
        }

        //! Performs an Y-Pauli gate
        //! \param qubit the qubit to perform the gate on.
        void y(int qubit) {
            maybe_single_error(qubit);
            baseSim.y(qubit);
        }

        //! Performs an Z-Pauli gate
        //! \param qubit the qubit to perform the gate on.
        void z(int qubit) {
            maybe_single_error(qubit);
            baseSim.z(qubit);
        }

        //! Performs an identity gate on a given qubit.
        //! \param qubit the qubit to perform the gate on.
        void idle(int qubit) {
            maybe_single_error(qubit);
            baseSim.idle(qubit);
        }

        //! Measures a qubit in the Z basis
        //! \param qubit
        //! \return The measurement outcome as a measurement_state;
        measurement_state measure(int qubit){
            maybe_measure_error_pre(qubit);
            auto result = baseSim.measure(qubit);
            maybe_measure_error_post(qubit);
            return result;
        }

        template<typename T>
        auto measure_all(T qubits){
            for (auto qubit: qubits){
                maybe_measure_error_pre(qubit);
            }
            auto retval = baseSim.measuere_all(qubits);
            for (auto qubit: qubits){
                maybe_measure_error_post(qubit);
            }
            return retval;
        }

        template<typename T>
        auto measure_reset_all(T qubits){
            for (auto qubit: qubits){
                maybe_measure_error_pre(qubit);
            }
            auto retval = baseSim.measuere_reset_all(qubits);
            for (auto qubit: qubits){
                maybe_measure_error_post(qubit);
            }
            return retval;
        }

        //! Measures a qubit in the Z basis and then resets it to single_pauli_state::ZERO.
        //! \param qubit
        //! \return The measurement outcome as a measurement_state;
        measurement_state measure_and_reset(int qubit) {
            maybe_measure_error_pre(qubit);
            auto result = baseSim.measure_and_reset(qubit);
            maybe_measure_error_post(qubit);
            return result;
        }

        //! Measures a qubit in the given pauli basis and then resets it to single_pauli_state::ZERO.
        //! \param qubit
        //! \param the basis to measure in
        //! \return The measurement outcome as a measurement_state;
        measurement_state measure_in_basis_and_reset(int qubit, single_pauli_basis b) {
            maybe_measure_error_pre(qubit);
            auto result = baseSim.measure_in_basis_and_reset(qubit, b);
            maybe_measure_error_post(qubit);
            return result;
        }

        //! Measures a qubit in the Z basis and then resets it to the state s;
        //! \param qubit
        //! \return The measurement outcome as a measurement_state;
        measurement_state measure_and_reset_to_state(int qubit, single_pauli_state s) {
            maybe_measure_error_pre(qubit);
            auto result = baseSim.measure_and_reset_to_state(qubit, s);
            maybe_measure_error_post(qubit);
            return result;
        }

        //! Measures a qubit in the given pauli basis and then resets it to single_pauli_state::ZERO.
        //! \param qubit
        //! \param the basis to measure in
        //! \return The measurement outcome as a measurement_state;
        measurement_state measure_in_basis_and_reset_to_state(int qubit, single_pauli_basis b, single_pauli_state s) {
            maybe_measure_error_pre(qubit);
            auto result = baseSim.measure_in_basis_and_reset_to_state(qubit, b, s);
            maybe_measure_error_post(qubit);
            return result;
        }

        //! Converts the state to a std::string that represents the tableu of the state.
        explicit operator std::string() {
            return std::string(baseSim);
        }

        //! compares equality using only the stabilizer part of the tableu, does not normalize in the comparison.
        //! in future we should change this to compare based on whether the tables represent the same states by gaussian
        //! elimination.
        inline bool operator==(FixedErrorSim <BaseSim> &other) {
            return baseSim == other.baseSim;
        }

        //! Compares if the stabiliziser parts are not equal.
        inline bool operator!=(FixedErrorSim <BaseSim> &other) {
            return !((*this) == other);
        }

        //! Helper function that performs a rotation of a given qubit from the ZERO state to the given state;
        //! This fuction does not add errors.
        //! \param qubit - The qubit to rotate
        //! \param state - The state to rotate to
        inline void rotate_qubit(int qubit, single_pauli_state state) {
            return baseSim.rotate_qubit(qubit, state);
        }

        //! Function that returns an error free copy of the base class, to enable clean compsition of Error correction simulations.
        error_free_sim &get_error_free() {
            return baseSim.get_error_free();
        }

        //! Function that returns an error free copy of the base class, to enable clean compsition of Error correction simulations.
        const error_free_sim &get_error_free() const {
            return baseSim.get_error_free();
        }
    };
}
#endif //STEANE_SIM_FIXED_ERROR_SIM_H
