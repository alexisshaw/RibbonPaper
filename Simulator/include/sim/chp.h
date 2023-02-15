//
// Created by Alexis Shaw on 25/05/2020.
//
// Based on the Python CHP Stabilizer Simulator Library,
// https://github.com/Strilanc/python-chp-stabilizer-simulator
//
// And the Algorithm described in "Improved Simulation of Stabilizer Circuits"
// by Scott Aaronson, and Daniel Gottesman "https://arxiv.org/abs/quant-ph/0406196"
//
// Licensed under the Apache License.
//

#ifndef _chp_simulator
#define _chp_simulator

#include <string>
#include <array>
#include <random>
#include <bitset>
#include <memory>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include "BaseTypes.h"
#include "prettyprint.h"
#include "magic_enum.h"

using namespace magic_enum::ostream_operators;

namespace CodeTools {
    //! A basic unoptimised implementation of the CHP algorithm proposed by Scott Aaronson
    //! The table is a 2d 2n+1 x 2n c-array of bools, with a separate vector for the stabilizer values.
    //! \tparam n An integer template parameter that contains the number of
    //! \tparam R The PRNG type used for measurements, by default sdd::mt19937_64.
    using namespace std::numbers;

    template<int n, typename R=std::mt19937_64>
    class ChpSimulator {
    public:
        static constexpr int N{n};
        static constexpr double epsilon = 1E-10;

        //! A base simulator class, operations on which occur error free
        using error_free_sim = ChpSimulator;

        //! This constructor takes in a shared_ptr to a random generator to allow the reuse of the random generator.
        explicit ChpSimulator(std::shared_ptr<R> randGenIn);

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

                //! Performs an X-Pauli gate
        //! \param qubit the qubit to perform the gate on.
        void x(int qubit, double alpha){
            if (std::abs(alpha) < epsilon){
                return;
            } else if (std::abs(alpha - pi) < epsilon){
                x(qubit);
                return;
            } else {
                throw std::runtime_error(std::string("impossible to perform rotation with arbitrary angle"));
            }
        };

        //! Performs an Y-Pauli gate
        //! \param qubit the qubit to perform the gate on.
        void y(int qubit, double alpha){
            if (std::abs(alpha) < epsilon){
                return;
            } else if (std::abs(alpha - pi) < epsilon){
                y(qubit);
                return;
            } else {
                throw std::runtime_error(std::string("impossible to perform rotation with arbitrary angle"));
            }
        };

        //! Performs an Z-Pauli gate
        //! \param qubit the qubit to perform the gate on.
        void z(int qubit, double alpha){
            if (std::abs(alpha) < epsilon){
                return;
            } else if (std::abs(alpha - pi) < epsilon){
                z(qubit);
                return;
            } else {
                throw std::runtime_error(std::string("impossible to perform rotation with arbitrary angle"));
            }
        };

        //! Performs an identity gate on a given qubit. In the base simulator this does nothing.
        //! \param qubit the qubit to perform the gate on.
        void idle(int qubit){};

        //! Measures a qubit in the Z basis.
        //! \param qubit
        //! \return The measurement outcome as a measurement_state;
        measurement_state measure(int qubit, bool rand = true, bool value = false);

        //! Measures a qubit in the Z basis and then resets it to single_pauli_state::ZERO.
        //! \param qubit
        //! \return The measurement outcome as a measurement_state;
        measurement_state measure_and_reset(int qubit, bool rand = true, bool value = false);

        //! Measures a qubit in the given pauli basis and then resets it to single_pauli_state::ZERO.
        //! \param qubit
        //! \param the basis to measure in
        //! \return The measurement outcome as a measurement_state;
        measurement_state measure_in_basis_and_reset(int qubit, single_pauli_basis b);

        //! Measures a qubit in the Z basis and then resets it to the state s;
        //! \param qubit
        //! \return The measurement outcome as a measurement_state;
        measurement_state measure_and_reset_to_state(int qubit, single_pauli_state s);

        //! Measures a qubit in the given pauli basis and then resets it to single_pauli_state::ZERO.
        //! \param qubit
        //! \param the basis to measure in
        //! \return The measurement outcome as a measurement_state;
        measurement_state measure_in_basis_and_reset_to_state(int qubit, single_pauli_basis b, single_pauli_state s);

        //! Converts the state to a std::string that represents the tableu of the state.
        explicit operator std::string();

        std::string row_as_string(int row);

        //! compares equality using only the stabilizer part of the tableu, does not normalize in the comparison.
        //! in future we should change this to compare based on whether the tables represent the same states by gaussian
        //! elimination.
        inline bool operator==(const ChpSimulator &other);

        //! Compares if the stabiliziser parts are not equal.
        inline bool operator!=(const ChpSimulator &other);

        //! Gets the inner product of two chp simulators, DUMMY implementation now.
        inline double getInnerProduct(const ChpSimulator& other){
            if((*this) == other){
                return 1.0;
            } else {
                return 0.0;
            }
        }

        //! A function used to help verify the inner portion of row_mult. Not actually called anywhere else.
        static int pauli_product_phase(bool x1, bool z1, bool x2, bool z2);

        //! Helper function that performs a rotation of a given qubit from the ZERO state to the given state;
        //! \param qubit - The qubit to rotate
        //! \param state - The state to rotate to
        void rotate_qubit(int qubit, single_pauli_state state);

        //! Function that returns an error free copy of the base class, to enable clean composition of Error correction simulations.
        error_free_sim& get_error_free(){
            return (*this);
        }

        //! Function that returns an error free copy of the base class, to enable clean composition of Error correction simulations.
        const error_free_sim& get_error_free() const{
            return (*this);
        }

        //! Function that applies an injected single pauli error, effectively just applying the pauli here, used to help implement
        //! error model adaptor template classes.
        void inject_single_pauli_error(int qubit, single_qubit_pauli p){
            switch (p) {
                case single_qubit_pauli::I:
                    break;
                case single_qubit_pauli::Z:
                    z(qubit);
                    break;
                case single_qubit_pauli::X:
                    x(qubit);
                    break;
                case single_qubit_pauli::Y:
                    y(qubit);
                    break;
            }
        }

        void normalize(void){
            gaussian();
            ut_gaussian();
            gaussian();
        }
        int gaussian(void);
        int ut_gaussian(void);
    protected:
        //! This contains the table of the pauli matrix, it is setup in the order such that table[qubit_no][stab_no]
        //! there are 4 quadrants,
        //! for 0 <= qubit_no < n, the row represents X part of (de)stabilisers for a given qubit.
        //! for n <= qubit_no < 2*n, the row represents Z part of (de)stabilisers for a given qubit.
        //! for 0 <= stab_no < n, the col represents a destabiliser of the state (as defined in the paper)
        //! for n <= stab_no < 2*n, the col represents a stabiliser of the state.
        std::array<std::bitset<2 * n>, 2 * n> table{};

        //! This contains the parities for the stabilizer matrix, it is setup in the order such that r[stab_no] is the parity
        //! for a given stabiliser. It is a bitset to help with fast arithmetic.
        //! for 0 <= stab_no < n, the represents a destabiliser of the state (as defined in the paper)
        //! for n <= stab_no < 2*n, the represents a stabiliser of the state.
        std::bitset<2 * n> r{};
        std::shared_ptr<R> rand_gen;

        //! Multiply two stabilisers in the table, given by their stab no.
        void row_mult(int i, int k);

        //! multiply two stabilisers, one given by references, and the other by reference into the table.
        //! return the value of the parity correction to be applied to the parity.
        bool row_mult(std::bitset<n> &x1, std::bitset<n> &z1, size_t k);

        std::string row_as_string(std::bitset<n> &x1, std::bitset<n> &z1, bool r);

        //! multiply two stabilisers, both given by references to bitsets describing their value.
        //! return the value of the parity correction to be applied to the parity.
        bool row_mult(std::bitset<n> &x1, std::bitset<n> &z1, std::bitset<n> &x2, std::bitset<n> &z2);

        //! swaps two different stabilisers in the table, given by their stab nos.
        void row_swap(int i, int k);

        //! Helper function for measuring a qubit that returns a random value.
        //! returns a measurement state.
        measurement_state measure_random(int qubit, int p_dash, bool rand=true, bool val = false);

        //! Helper functon for measuring a qubit that returns a deterministic value.
        //! returns a measurement state.
        measurement_state measure_deterministic(int qubit);

        //! Helper function used in the gauusian elimination routine.
        void do_elimination(std::bitset<2*n> &b, int &i);
        void do_reverse_elimination(std::bitset<2*n> &b, int &i);

        bool check_span(std::bitset<n> &x, std::bitset<n> &z, bool r_val, std::array<size_t, n> &Indexes);

        public:
        template<typename T>
        boost::container::small_vector<measurement_state, 20> measure_all(T qubits){
            boost::container::small_vector<measurement_state, 20> retval{};
            for (auto qubit_no: qubits){
                retval.push_back(measure(qubit_no));
            }
            return retval;
        }

        template<typename T>
        boost::container::small_vector<measurement_state, 20> measure_reset_all(T qubits){
            boost::container::small_vector<measurement_state, 20> retval{};
            for (auto qubit_no: qubits){
                retval.push_back(measure_and_reset(qubit_no));
            }
            return retval;
        }
        
        template<typename T, typename U>
        double measure_all_set(T qubits, U required){
            boost::container::small_vector<measurement_state, 20> measurement_results{};
            for (size_t i = 0; i < qubits.size(); ++i){
                bool value = (required[i] == measurement_state::DETERMINISTIC_TRUE) || (required[i] == measurement_state::RANDOM_TRUE);
                measurement_results.push_back(measure(qubits[i], false, value));
                // std::cout << required[i] << ' ' <<measurement_results[i] << std::endl;
            }
            

            double p = 1;
            for (size_t i = 0; i < qubits.size(); ++i){
                if(!measurement_state_random(measurement_results[i]) && 
                   (measurement_state_value(measurement_results[i]) != measurement_state_value(required[i]))){
                        p = 0.0;
                }
                if(measurement_state_random(measurement_results[i])){
                    if (measurement_results[i] != required[i]){
                        std::cerr << "incorrect value for pre-determined random measurement error" << std::endl; 
                        abort();
                    }
                    p = p * 0.5;
                }
            }
            return p;
        }

        template<typename T, typename U>
        double measure_reset_all_set(T qubits, U required){
            boost::container::small_vector<measurement_state, 20> measurement_results{};
            for (size_t i = 0; i < qubits.size(); ++i){
                bool value = (required[i] == measurement_state::DETERMINISTIC_TRUE) || (required[i] == measurement_state::RANDOM_TRUE);
                measurement_results.push_back(measure_and_reset(qubits[i], false, value));
                // std::cout << required[i] << ' ' <<measurement_results[i] << std::endl;
            }
            

            double p = 1;
            for (size_t i = 0; i < qubits.size(); ++i){
                if(!measurement_state_random(measurement_results[i]) && 
                   (measurement_state_value(measurement_results[i]) != measurement_state_value(required[i]))){
                        p = 0.0;
                }
                if(measurement_state_random(measurement_results[i])){
                    if (measurement_results[i] != required[i]){
                        std::cerr << "incorrect value for pre-determined random measurement error" << std::endl; 
                        abort();
                    }
                    p = p * 0.5;
                }
            }
            return p;
        }
    };

    template<int n, typename R>
    inline ChpSimulator<n, R>::ChpSimulator(std::shared_ptr<R> randGenIn):rand_gen{randGenIn}{
        // Setup Table
        for (int i = 0; i < 2 * n; i++) {
            table[i][i] = true;
        }
    }

    template<int n, typename R>
    inline bool ChpSimulator<n, R>::operator==(const ChpSimulator &other) {
        bool initial_test = true;
        [&]{
        for (int i = 0; i != 2 * n; i++)
            for (int j = 0; j < n; j++)
                if (table[i][j + n] != other.table[i][j + n]){
                    initial_test = false;
                    return;
                }
        for(int i = 0; i != 2*n; i++)
            if(r[i+n] != other.r[i + n]){
                initial_test = false;
                return;
            }
        }();
        if(initial_test){
            return true;
        }
        ChpSimulator temp_this{*this};

        temp_this.gaussian();
        std::array<size_t, n> indexes{};
        // std::cout << std::string(temp_this) << std::endl;

        for (int i = 0; i < n; i++){
            indexes[i] = std::find_if(begin(temp_this.table),end(temp_this.table),
                                   [&](auto &j){return j[i + n];}) - begin(temp_this.table);
        }
        // std::cout << indexes << std::endl;

        for (int i = 0; i < n ; i++){
            std::bitset<n> x{}, z{};
            for (int j = 0; j < n; j++) {
                x[j] = other.table[j][i+n];
                z[j] = other.table[j + n][i+n];
            }
            if (!temp_this.check_span(x, z, other.r[i+n] , indexes))
                return false;
        }

        return true;
    }

    template<int n, typename R>
    inline bool ChpSimulator<n, R>::operator!=(const ChpSimulator &other) {
        return !((*this) == other);
    }

    template<int n, typename R>
    inline void ChpSimulator<n, R>::h(int qubit) {
        auto &x = table[qubit];
        auto &z = table[qubit + n];
        auto &parity = r;

        parity = parity ^ (x & z);

        // swap elements.
        std::swap(z, x);
    }

    template<int n, typename R>
    inline void ChpSimulator<n, R>::s(int qubit) {
        auto &x = table[qubit];
        auto &z = table[qubit + n];
        auto &parity = r;

        parity = parity ^ (x & z);
        z = z ^ x;
    }

    template<int n, typename R>
    inline void ChpSimulator<n, R>::sdg(int qubit) {
        auto &x = table[qubit];
        auto &z = table[qubit + n];
        auto &parity = r;

        parity = parity ^ (x & (~z));
        z = z ^ x;
    }


    template<int n, typename R>
    inline void ChpSimulator<n, R>::x(int qubit) {
        auto &z = table[qubit + n];
        auto &parity = r;

        parity = parity ^ z;
    }

    template<int n, typename R>
    inline void ChpSimulator<n, R>::y(int qubit) {
        auto &x = table[qubit];
        auto &z = table[qubit + n];
        auto &parity = r;

        parity = parity ^ x ^ z;
    }

    template<int n, typename R>
    inline void ChpSimulator<n, R>::z(int qubit) {
        auto &x = table[qubit];
        auto &parity = r;

        parity = parity ^ x;
    }

    template<int n, typename R>
    inline void ChpSimulator<n, R>::cx(int control, int target) {
        auto &x_c = table[control];
        auto &z_c = table[control + n];
        auto &x_t = table[target];
        auto &z_t = table[target + n];
        auto &parity = r;

        parity = parity ^ (x_c & z_t & ~(x_t ^ z_c));
        x_t = x_t ^ x_c;
        z_c = z_c ^ z_t;
    }

    template<int n, typename R>
    inline void ChpSimulator<n, R>::cy(int control, int target) {
        auto &x_c = table[control];
        auto &z_c = table[control + n];
        auto &x_t = table[target];
        auto &z_t = table[target + n];
        auto &parity = r;

        parity = parity ^ (x_c & ((~x_t & z_c & z_t) | (x_t & ~z_c & ~z_t)));
        z_c = z_c ^ x_t ^ z_t;
        z_t = z_t ^ x_c;
        x_t = x_t ^ x_c;
    }

    template<int n, typename R>
    inline void ChpSimulator<n, R>::cz(int control, int target) {
        auto &x_c = table[control];
        auto &z_c = table[control + n];
        auto &x_t = table[target];
        auto &z_t = table[target + n];
        auto &parity = r;

        parity = parity ^ (x_c & x_t & (z_t ^ z_c));
        z_t = z_t ^ x_c;
        z_c = z_c ^ x_t;
    }

    template<typename T>
    inline int get_next(const T &array, int start) {
        for (; start < (static_cast<int>(array.size())); ++start) {
            if (array[start])
                return start;
        }
        return start;
    }

    template<int n, typename R>
    inline measurement_state ChpSimulator<n, R>::measure(int qubit, bool rand, bool val) {
        // find_first.
        int p = get_next(table[qubit], n);
        if (p < 2 * n)
            return measure_random(qubit, p, rand, val);
        return measure_deterministic(qubit);
    }

    template<int n, typename R>
    inline measurement_state ChpSimulator<n, R>::measure_random(int qubit, int p, bool rand, bool val) {
        int p_dash = p - n;
        r[p_dash] = r[p];
        for (int i = 0; i < (2 * n); i++) {
            table[i][p_dash] = bool(table[i][p]);
            table[i][p] = false;
        }

        table[qubit + n][p] = true;
        bool rand_val = val;
        if(rand){
            rand_val = ((*rand_gen)() > ((rand_gen->max() - rand_gen->min()) / 2));
        }
        r[p] = rand_val;

        for (int i = 0; i < (2 * n); i++) {
            if ((i != p_dash) & (i != p) & (table[qubit][i])) {
                row_mult(i, p_dash);
            }
        }

        return r[p] ? measurement_state::RANDOM_TRUE : measurement_state::RANDOM_FALSE;
    }

    template<int n, typename R>
    inline measurement_state ChpSimulator<n, R>::measure_deterministic(int qubit) {
        // Temporary Storage of an additional row for a deterministic measurement.
        std::bitset<n> x1{}, z1{};
        bool r_2n = false;
        for (int i = 0; i < n; i++) {
            if (table[qubit][i]) {
                bool p = row_mult(x1, z1, i + n);
                r_2n = bool(r[i + n] ^ p ^ r_2n);
            }
        }

        return r_2n ? measurement_state::DETERMINISTIC_TRUE : measurement_state::DETERMINISTIC_FALSE;
    }

    template<int n, typename R>
    inline measurement_state ChpSimulator<n, R>::measure_and_reset(int qubit, bool rand, bool value) {
        auto measurement = measure(qubit, rand, value);
        if ((measurement == measurement_state::DETERMINISTIC_TRUE) | (measurement == measurement_state::RANDOM_TRUE)) {
            x(qubit);
        }
        return measurement;
    }

    template<int n, typename R>
    inline measurement_state ChpSimulator<n, R>::measure_and_reset_to_state(int qubit, single_pauli_state state) {
        auto measurement = measure(qubit);
        if ((measurement == measurement_state::DETERMINISTIC_TRUE) | (measurement == measurement_state::RANDOM_TRUE)) {
            x(qubit);
        }
        switch (state) {
            case single_pauli_state::ZERO:
                break;
            case single_pauli_state::PLUS:
                h(qubit);
                break;
            case single_pauli_state::I_PLUS:
                h(qubit);
                s(qubit);
                break;
            case single_pauli_state::ONE:
                x(qubit);
                break;
            case single_pauli_state::MINUS:
                x(qubit);
                h(qubit);
                break;
            case single_pauli_state::I_MINUS:
                x(qubit);
                h(qubit);
                s(qubit);
            default:
                break;
        }
        return measurement;
    }

    template<int n, typename R>
    inline measurement_state ChpSimulator<n, R>::measure_in_basis_and_reset_to_state(int qubit, single_pauli_basis basis, single_pauli_state state) {
        switch(basis){
            case single_pauli_basis::Z:
                break;
            case single_pauli_basis::X:
                h(qubit);
                break;
            case single_pauli_basis::Y:
                sdg(qubit);
                h(qubit);
                break;
            default:
                break;
        }

        auto measurement = measure(qubit);
        if ((measurement == measurement_state::DETERMINISTIC_TRUE) | (measurement == measurement_state::RANDOM_TRUE)) {
            x(qubit);
        }
        switch (state) {
            case single_pauli_state::ZERO:
                break;
            case single_pauli_state::PLUS:
                h(qubit);
                break;
            case single_pauli_state::I_PLUS:
                h(qubit);
                s(qubit);
                break;
            case single_pauli_state::ONE:
                x(qubit);
                break;
            case single_pauli_state::MINUS:
                x(qubit);
                h(qubit);
                break;
            case single_pauli_state::I_MINUS:
                x(qubit);
                h(qubit);
                s(qubit);
            default:
                break;
        }
        return measurement;
    }

    template<int n, typename R>
    inline measurement_state ChpSimulator<n, R>::measure_in_basis_and_reset(int qubit, single_pauli_basis basis) {
        switch(basis){
            case single_pauli_basis::Z:
                break;
            case single_pauli_basis::X:
                h(qubit);
                break;
            case single_pauli_basis::Y:
                sdg(qubit);
                h(qubit);
                break;
            default:
                break;
        }

        auto measurement = measure(qubit);
        if ((measurement == measurement_state::DETERMINISTIC_TRUE) | (measurement == measurement_state::RANDOM_TRUE)) {
            x(qubit);
        }
        return measurement;
    }

    template<int n, typename R>
    inline int
    ChpSimulator<n,R>::gaussian(){
        int i = n;
        for (int j = 0; j < n; j++){
            do_elimination(table[j], i);
        }

        int g = i - n;

        for (int j = 0; j < n; j++){
            do_elimination(table[j+n], i);
        }
        return g;
    }

    template<int n, typename R>
    void ChpSimulator<n,R>::do_elimination(std::bitset<2*n> &b, int &i){
        int k;
        // Find a generator containing set bit after i
        for (k = i; k < 2*n; k++) 
                    if (b[k]) break;
        if (k < 2*n){
            row_swap(i, k);
            row_swap(i-n, k-n);
            for (int k2 = i+1; k2 < 2*n; k2++)
                if (b[k2]){
                    row_mult(k2, i);
                    row_mult(i-n, k2-n);
                }
            i++;
        }
    }

    template<int n, typename R>
    inline int
    ChpSimulator<n,R>::ut_gaussian(){
        int i = 2*n - 1;
        for (int j = (2 * n-1); j >= 0; j--){
            do_reverse_elimination(table[j], i);
        }
        return 0;
    }

    template<int n, typename R>
    void ChpSimulator<n,R>::do_reverse_elimination(std::bitset<2*n> &b, int &i){
        int k;
        // Find a generator containing set bit after i
        for (k = i; k >= n; k--) 
                    if (b[k]) break;
        // std::cout << k << " " << b << " " << i << std::endl;
        if (k >= n){
            row_swap(i, k);
            row_swap(i-n, k-n);
            for (int k2 = i-1; k2 >= n; k2--)
                if (b[k2]){
                    row_mult(k2, i);
                    row_mult(i-n, k2-n);
                }
            i--;
        }
    }

    template<int n, typename R>
    bool ChpSimulator<n,R>::check_span(std::bitset<n> &x, std::bitset<n> &z, bool r_val, std::array<size_t, n> &indexes){
        // Assumes that this row has been put into upper-triangular form 
        // Indexes contains the index of the first non-zero bit in the associated syndrome
        
        // std::cout << "====================\nchecking syndrome:\n"<< row_as_string(x,z,r_val) <<"====================\n"<< std::flush;

        for (size_t i = 0; i < static_cast<size_t>(n); ++i){
            if(x[i]){
                // std::cout << "looking for stab with X/Y value in collumn " << i << std::endl;
                auto result = std::lower_bound(begin(indexes), end(indexes), i);
                if((*result) != i){
                    // std::cout << "not found" << std::endl;
                    return false; //There is no element that can eliminate the element.
                }
                size_t stab_id = result - begin(indexes) + n;
                // std::cout << "multiplying with stabiliser no: " << stab_id << "\n" << row_as_string(stab_id) << std::flush;
                r_val = r_val ^ r[stab_id] ^ row_mult(x,z, stab_id);
                // std::cout << "--------\n" << row_as_string(x,z,r_val) << std::endl;
            }
        }
        for (size_t i = 0; i < static_cast<size_t>(n); ++i){
            if(z[i]){
                // std::cout << "looking for stab with Z value in collumn " << i << " (" << i+n << ")" << std::endl;
                auto result = std::lower_bound(begin(indexes), end(indexes), i + n);                
                if((*result) != (i+n)){
                    // std::cout << "not found" << std::endl;
                    return false; //There is no element that can eliminate the element.
                }
                size_t stab_id = result - begin(indexes) + n;
                // std::cout << "multiplying with stabiliser no: " << stab_id << "\n" << row_as_string(stab_id) << std::flush;
                r_val = r_val ^ r[stab_id] ^ row_mult(x,z, stab_id);
                // std::cout << "--------\n" << row_as_string(x,z,r_val) << std::endl;
            }
        }
        // std::cout << std::endl;
        return !r_val;
    }


    template<int n, typename R>
    inline bool
    ChpSimulator<n, R>::row_mult(std::bitset<n> &x1, std::bitset<n> &z1, std::bitset<n> &x2, std::bitset<n> &z2) {
        size_t pauli_phases = 0u;

        std::bitset<n> neg{(~x1 & z1 & x2 & ~z2) | (x1 & z2 & (z1 ^ x2))};
        std::bitset<n> pos{(x1 & z2 & ~z1 & ~x2) | (z1 & x2 & (x1 ^ z2))};

        pauli_phases = neg.count() - pos.count();
        x1 ^= x2;
        z1 ^= z2;
        return ((pauli_phases & 2u) >> 1u);
    }

    template<int n, typename R>
    inline bool ChpSimulator<n, R>::row_mult(std::bitset<n> &x1, std::bitset<n> &z1, size_t k) {
        std::bitset<n> x2{}, z2{};

        for (int j = 0; j < n; j++) {
            x2[j] = table[j][k];
            z2[j] = table[j + n][k];
        }

        return row_mult(x1, z1, x2, z2);
    }

    template<int n, typename R>
    inline void ChpSimulator<n, R>::row_mult(int i, int k) {
        std::bitset<n> x1{}, z1{}, x2{}, z2{};

        for (int j = 0; j < n; j++) {
            x1[j] = table[j][i];
            x2[j] = table[j][k];
            z1[j] = table[j + n][i];
            z2[j] = table[j + n][k];
        }

        bool p = row_mult(x1, z1, x2, z2);
        r[i] = bool(r[i] ^ r[k] ^ p);

        for (int j = 0; j < n; ++j) {
            table[j][i] = x1[j];
        }
        for (int j = 0; j < n; ++j) {
            table[j + n][i] = z1[j];
        }
    }

    template<int n, typename R>
    inline void ChpSimulator<n, R>::row_swap(int i, int k) {
        for (int j = 0; j < n*2; j++) {
            bool t  = table[j][i];
            table[j][i] = table[j][k];
            table[j][k] = t;
        }

        bool t = r[i];
        r[i] = r[k];
        r[k]= t;
    }



    template<int n, typename R>
    inline int ChpSimulator<n, R>::pauli_product_phase(bool x1, bool z1, bool x2, bool z2) {
        //Determines the power of i in the product of two Paulis.
        //For example, X*Y = iZ and so this method would return +1 for X and Y.
        //The input Paulis are encoded into the following form:
        //x z | Pauli
        //----+-------
        //0 0 | I
        //1 0 | X
        //1 1 | Y
        //0 1 | Z
        //
        //Y gate:
        //No phase for YI = Y
        //-1 phase for YX = -iZ
        //No phase for YY = I
        //+1 phase for YZ = +iX
        //
        //X gate.
        //No phase for XI = X
        //No phase for XX = I
        //+1 phase for XY = iZ
        //-1 phase for XZ = -iY
        //
        //Z gate.
        //No phase for ZI = Z
        //+1 phase for ZX = -iY
        //-1 phase for ZY = iX
        //No phase for ZZ = I
        //
        //This is then done as a pauli basis

        bool neg = (!x1 && z1 && x2 && !z2) | (x1 && z2 && (z1 ^ x2));
        bool pos = (x1 && z2 && !z1 && !x2) | (z1 && x2 && (x1 ^ z2));

        return int(neg) - int(pos);
    }

    template<int n, typename R>
    std::string ChpSimulator<n, R>::row_as_string(int row){
        std::string ret{};
        const char outputs[4] = {'.', 'X', 'Z', 'Y'};
        ret.push_back(r[row] ? '-' : '+');
        for (int col = 0; col < n; col++) {
            int k = 2 * table[col + n][row] + table[col][row];
            ret.push_back(outputs[k]);
        }
        ret.push_back('\n');
        return ret;
    }

    template<int n, typename R>
    std::string ChpSimulator<n, R>::row_as_string(std::bitset<n> &x, std::bitset<n> &z, bool r_in){
        std::string ret{};
        const char outputs[4] = {'.', 'X', 'Z', 'Y'};
        ret.push_back(r_in ? '-' : '+');
        for (int col = 0; col < n; col++) {
            int k = 2 * z[col] + x[col];
            ret.push_back(outputs[k]);
        }
        ret.push_back('\n');
        return ret;
    }

    template<int n, typename R>
    inline ChpSimulator<n, R>::operator std::string() {

        auto ret = std::string();
        for (int row = 0; row < 2 * n; row++) {
            if (row == n) {
                for (int col = 0; col < (n + 1); col++)
                    ret.push_back('-');
                ret.push_back('\n');
            }

            ret += row_as_string(row);
        }
        return ret;
    }

    template<int n, typename R>
    inline void ChpSimulator<n, R>::rotate_qubit(const int qubit, single_pauli_state state) {
        // ZERO, PLUS, I_PLUS, ONE, MINUS, I_MINUS
        switch (state) {
            case single_pauli_state::ZERO:
                return;
            case single_pauli_state::PLUS:
                h(qubit);
                return;
            case single_pauli_state::I_PLUS:
                h(qubit);
                s(qubit);
                return;
            case single_pauli_state::ONE:
                x(qubit);
                return;
            case single_pauli_state::MINUS:
                x(qubit);
                h(qubit);
                return;
            case single_pauli_state::I_MINUS:
                x(qubit);
                h(qubit);
                s(qubit);
            default:
                return;
        }
    }
}

#endif