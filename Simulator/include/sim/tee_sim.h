//
// Created by 13114347 on 22/07/2020.
//

#ifndef STEANE_SIM_TEE_H
#define STEANE_SIM_TEE_H

#include <string>
#include <random>
#include "BaseTypes.h"

namespace CodeTools {
    template<typename BaseSim, typename AltSim>
    class TeeSim {
    public:
        static constexpr int N{std::min(BaseSim::N, AltSim::N)};

        BaseSim baseSim;
        AltSim altSim;

        //! A base simulator class, operations on which occur error free
        using error_free_sim = typename BaseSim::error_free_sim;

        //! Creates a n
        TeeSim(const BaseSim &b, const AltSim &a) : baseSim{b}, altSim{a} {};

        void inject_single_pauli_error(int qubit, single_qubit_pauli p) {
            baseSim.inject_single_pauli_error(qubit, p);
            altSim.inject_single_pauli_error(qubit, p);
        };

        //! Enables an error model sim to begin a moment;
        //! Does nothing here
        void moment(){
            altSim.moment();
            baseSim.moment();
        };

        //! Performs a controlled-X gate, (CNOT gate).
        //! \param control The control qubit id.
        //! \param target The target qubit id.
        void cx(int control, int target) {
            altSim.cx(control,target);
            baseSim.cx(control, target);
        }

        //! Performs a controlled-Y gate
        //! \param control The control qubit id.
        //! \param target The target qubit id.
        void cy(int control, int target) {
            altSim.cy(control, target);
            baseSim.cy(control, target);
        }

        //! Performs a controlled-Z gate
        //! \param control The control qubit id.
        //! \param target The target qubit id.
        void cz(int control, int target) {
            altSim.cz(control, target);
            baseSim.cz(control, target);
        }

        //! Performs a Hadamard Gate
        //! \param qubit the qubit to perform the gate on.
        void h(int qubit) {
            altSim.h(qubit);
            baseSim.h(qubit);
        }

        //! Performs a Phase Gate
        //! \param qubit the qubit to perform the gate on.
        void s(int qubit) {
            altSim.s(qubit);
            baseSim.s(qubit);
        }

        //! Performs an inverse Phase Gate
        //! \param qubit the qubit to perform the gate on.
        void sdg(int qubit) {
            altSim.sdg(qubit);
            baseSim.sdg(qubit);
        }

        //! Performs an X-Pauli gate
        //! \param qubit the qubit to perform the gate on.
        void x(int qubit) {
            altSim.x(qubit);
            baseSim.x(qubit);
        }

        //! Performs an Y-Pauli gate
        //! \param qubit the qubit to perform the gate on.
        void y(int qubit) {
            altSim.y(qubit);
            baseSim.y(qubit);
        }

        //! Performs an Z-Pauli gate
        //! \param qubit the qubit to perform the gate on.
        void z(int qubit) {
            altSim.z(qubit);
            baseSim.z(qubit);
        }

        //! Performs an identity gate on a given qubit.
        //! \param qubit the qubit to perform the gate on.
        void idle(int qubit) {
            altSim.idle(qubit);
            baseSim.idle(qubit);
        }

        //! Measures a qubit in the Z basis
        //! \param qubit
        //! \return The measurement outcome as a measurement_state;
        measurement_state measure(int qubit){
            altSim.measure(qubit);
            return baseSim.measure(qubit);
        }

        template<typename T>
        auto measure_all(T qubits){
            altSim.measure_all(qubits);
            return baseSim.measure_all(qubits);
        }

        //! Measures a qubit in the Z basis and then resets it to single_pauli_state::ZERO.
        //! \param qubit
        //! \return The measurement outcome as a measurement_state;
        measurement_state measure_and_reset(int qubit) {
            altSim.measure_and_reset(qubit);
            return baseSim.measure_and_reset(qubit);
        }

        //! Measures a qubit in the given pauli basis and then resets it to single_pauli_state::ZERO.
        //! \param qubit
        //! \param the basis to measure in
        //! \return The measurement outcome as a measurement_state;
        measurement_state measure_in_basis_and_reset(int qubit, single_pauli_basis b) {
            altSim.measure_in_basis_and_reset(qubit, b);
            return baseSim.measure_in_basis_and_reset(qubit, b);
        }

        //! Measures a qubit in the Z basis and then resets it to the state s;
        //! \param qubit
        //! \return The measurement outcome as a measurement_state;
        measurement_state measure_and_reset_to_state(int qubit, single_pauli_state s) {
            altSim.measure_and_reset_to_state(qubit, s);
            return baseSim.measure_and_reset_to_state(qubit, s);
        }

        //! Measures a qubit in the given pauli basis and then resets it to single_pauli_state::ZERO.
        //! \param qubit
        //! \param the basis to measure in
        //! \return The measurement outcome as a measurement_state;
        measurement_state measure_in_basis_and_reset_to_state(int qubit, single_pauli_basis b, single_pauli_state s) {
            altSim.measure_in_basis_and_reset_to_state(qubit, b, s);
            return baseSim.measure_in_basis_and_reset_to_state(qubit, b, s);
        }

        //! Converts the state to a std::string that represents the tableu of the state.
        explicit operator std::string() {
            return std::string(altSim) + "\n" + std::string(baseSim);
        }

        //! compares equality using only the stabilizer part of the tableu, does not normalize in the comparison.
        //! in future we should change this to compare based on whether the tables represent the same states by gaussian
        //! elimination.
        inline bool operator==(BalancedErrorSim <BaseSim> &other) {
            return baseSim == other.baseSim;
        }

        //! Compares if the stabiliziser parts are not equal.
        inline bool operator!=(BalancedErrorSim <BaseSim> &other) {
            return !((*this) == other);
        }

        //! Helper function that performs a rotation of a given qubit from the ZERO state to the given state;
        //! This fuction does not add errors.
        //! \param qubit - The qubit to rotate
        //! \param state - The state to rotate to
        inline void rotate_qubit(int qubit, single_pauli_state state) {
            altSim.rotate_quibit(qubit, state);
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

#endif //STEANE_SIM_TEE_H
