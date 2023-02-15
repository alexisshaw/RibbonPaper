//
// Created by Alexis Shaw on 24/06/2020.
//
//
// Licensed under the Apache License.
//

#ifndef STEANE_SIM_ERRORTRACKER_H
#define STEANE_SIM_ERRORTRACKER_H

#include <string>
#include <array>
#include <vector>
#include "BaseTypes.h"
namespace CodeTools{
    //! An error propagation tracker, to be used to determine appropriate syndromes for extraction, does NOT
    //! include any provision for accounting accumulated phase in injected errors.
    template<int n>
    class ErrorTracker {
    public:
        static constexpr int N{n};

        //! injects an error that is to be tracked.
        void inject_single_pauli_error(int qubit, single_qubit_pauli p);

        //! Enables an error model sim to begin a moment;
        //! Does nothing here
        void moment(){};

        //! Performs a controlled-X gate, (CNOT gate).
        //! \param control The control qubit id.
        //! \param target The target qubit id.
        void cx(int control, int target);

        //! Performs a controlled-Y gate
        //! \param control The control qubit id.
        //! \param target The target qubit id.
        void cy(int control, int target);

        //! Performs a controlled-Z gate
        //! \param control The control qubit id.
        //! \param target The target qubit id.
        void cz(int control, int target);

        //! Performs a Hadamard Gate
        //! \param qubit the qubit to perform the gate on.
        void h(int qubit);

        //! Performs a Phase Gate
        //! \param qubit the qubit to perform the gate on.
        void s(int qubit);

        //! Performs an inverse Phase Gate
        //! \param qubit the qubit to perform the gate on.
        void sdg(int qubit);

        //! Performs an X-Pauli gate
        //! \param qubit the qubit to perform the gate on.
        void x(int qubit);

        //! Performs an Y-Pauli gate
        //! \param qubit the qubit to perform the gate on.
        void y(int qubit);

        //! Performs an Z-Pauli gate
        //! \param qubit the qubit to perform the gate on.
        void z(int qubit);

        //! Performs an identity gate on a given qubit.
        //! \param qubit the qubit to perform the gate on.
        void idle(int qubit);

        //! An initialization on a qubit wipes away the error that was on that measurement, this facilitates that.
        void measure_and_reset(int qubit);

        //! Alias for consistance with SteaneSim
        void measure_in_basis_and_reset(int qubit, single_pauli_basis b);

        //! Alias for consistance with SteaneSim
        void measure_and_reset_to_state(int qubit, single_pauli_state s);

        //! Alias for consistance with SteaneSim
        void measure_in_basis_and_reset_to_state(int qubit, single_pauli_basis b, single_pauli_state s);

        //! Converts the state to a std::string that represents the error
        explicit operator std::string();

        error_correction_list getCorrection() const;
    protected:
        //! This contains the table of the pauli matrix, it is setup in the order such that table[qubit_no]
        //! has two halves.
        //! for 0 <= qubit_no < n, the row represents X part of (de)stabilisers for a given qubit.
        //! for n <= qubit_no < 2*n, the row represents Z part of (de)stabilisers for a given qubit.
        std::array<bool, n*2> table{};
        std::array<bool, n> dirty{};
    };

    template<int n>
    inline void ErrorTracker<n>::inject_single_pauli_error(int qubit, single_qubit_pauli p){
        auto &x = table[qubit];
        auto &z = table[qubit + n];
        switch(p){
            case single_qubit_pauli::I:
                break;
            case single_qubit_pauli::X:
                x = !x;
                break;
            case single_qubit_pauli::Y:
                x = !x;
                z = !z;
                break;
            case single_qubit_pauli::Z:
                z = !z;
                break;
        }
    }

    template<int n>
    inline void ErrorTracker<n>::h(int qubit) {
        auto &x = table[qubit];
        auto &z = table[qubit + n];
//        auto &parity = r;
//
//        parity = parity ^ (x & z);

        // swap elements.
        std::swap(z, x);
    }

    template<int n>
    inline void ErrorTracker<n>::s(int qubit) {
        auto &x = table[qubit];
        auto &z = table[qubit + n];
//        auto &parity = r;
//
//        parity = parity ^ (x & z);
        z = z ^ x;
    }

    template<int n>
    inline void ErrorTracker<n>::sdg(int qubit) {
        auto &x = table[qubit];
        auto &z = table[qubit + n];
//        auto &parity = r;
//
//        parity = parity ^ (x & (~z));
        z = z ^ x;
    }


    template<int n>
    inline void ErrorTracker<n>::x(int qubit) {
        auto &z = table[qubit + n];
//        auto &parity = r;
//
//        parity = parity ^ z;
    }

    template<int n>
    inline void ErrorTracker<n>::y(int qubit) {
        auto &x = table[qubit];
        auto &z = table[qubit + n];
//        auto &parity = r;
//
//        parity = parity ^ x ^ z;
    }

    template<int n>
    inline void ErrorTracker<n>::z(int qubit) {
        auto &x = table[qubit];
//        auto &parity = r;
//
//        parity = parity ^ x;
    }

    template<int n>
    inline void ErrorTracker<n>::cx(int control, int target) {
        auto &x_c = table[control];
        auto &z_c = table[control + n];
        auto &x_t = table[target];
        auto &z_t = table[target + n];
//        auto &parity = r;
//
//        parity = parity ^ (x_c & z_t & ~(x_t ^ z_c));
        x_t = x_t ^ x_c;
        z_c = z_c ^ z_t;
    }

    template<int n>
    inline void ErrorTracker<n>::cy(int control, int target) {
        auto &x_c = table[control];
        auto &z_c = table[control + n];
        auto &x_t = table[target];
        auto &z_t = table[target + n];
//        auto &parity = r;
//
//        parity = parity ^ (x_c & ((~x_t & z_c & z_t) | (x_t & ~z_c & ~z_t)));
        z_c = z_c ^ x_t ^ z_t;
        z_t = z_t ^ x_c;
        x_t = x_t ^ x_c;
    }

    template<int n>
    inline void ErrorTracker<n>::cz(int control, int target) {
        auto &x_c = table[control];
        auto &z_c = table[control + n];
        auto &x_t = table[target];
        auto &z_t = table[target + n];
//        auto &parity = r;
//
//        parity = parity ^ (x_c & x_t & (z_t ^ z_c));
        z_t = z_t ^ x_c;
        z_c = z_c ^ x_t;
    }

    template<int n>
    inline void ErrorTracker<n>::measure_and_reset(int qubit) {
        table[qubit] = false;
        table[qubit + n] = false;
    }

    template<int n>
    inline void ErrorTracker<n>::measure_and_reset_to_state(int qubit, single_pauli_state state) {
        measure_and_reset(qubit);
    }

    template<int n>
    inline void ErrorTracker<n>::measure_in_basis_and_reset_to_state(int qubit, single_pauli_basis basis, single_pauli_state state) {
        measure_and_reset(qubit);
    }

    template<int n>
    inline void ErrorTracker<n>::measure_in_basis_and_reset(int qubit, single_pauli_basis basis) {
        measure_and_reset(qubit);
    }

    template<int n>
    inline ErrorTracker<n>::operator std::string() {
        constexpr char outputs[4] = {'.', 'X', 'Z', 'Y'};
        auto ret = std::string();
        for (int col = 0; col < n; col++) {
            int k = 2 * table[col + n] + table[col];
            ret.push_back(outputs[k]);
        }
        return ret;
    }

    template<int n>
    inline error_correction_list ErrorTracker<n>::getCorrection() const{
        error_correction_list l{};
        for (int qubit = 0; qubit < n; qubit++){
            unsigned int k = 2*table[qubit + n] + table[qubit];
            switch ( k & 3u) {
                case 0:
                    // If the operator is an identity, then there is no operator to push back.
                    break;
                case 1:
                    l.push_back(std::make_pair(qubit, single_qubit_pauli::X));
                    break;
                case 2:
                    l.push_back(std::make_pair(qubit, single_qubit_pauli::Z));
                    break;
                case 3:
                    l.push_back(std::make_pair(qubit, single_qubit_pauli::Y));
                    break;
            }
        }
        return l;
    }

    template<int n>
    void ErrorTracker<n>::idle(int qubit) {
    }

}

#endif
