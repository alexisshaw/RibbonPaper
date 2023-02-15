
#ifndef _state_simulator
#define _state_simulator

#include <random>
#include <utility>
#include <memory>
#include "BaseTypes.h"
#include "sim/sim.h"
#include "boost/container/small_vector.hpp"

namespace CodeTools {
    //! A basic unoptimised implementation of the CHP algorithm proposed by Scott Aaronson
    //! The table is a 2d 2n+1 x 2n c-array of bools, with a separate vector for the stabilizer values.
    //! \tparam n An integer template parameter that contains the number of
    //! \tparam R The PRNG type used for measurements, by default sdd::mt19937_64.
    template<int n, typename R=std::mt19937_64>
    class StateSimulator {
    public:
        ds_Register reg;
        //int bit[1];

        static constexpr double PI_2 = 1.57079632679489661923;
    public:
        static constexpr int N{n};
        int enable_errors = 1;

        //! A base simulator class, operations on which occur error free
        using error_free_sim = StateSimulator;

//        static void initialise(){
//            ds_initialize_simulator();
//        }

        //! This constructor takes in a shared_ptr to a random generator to allow the reuse of the random generator.
        explicit StateSimulator(std::shared_ptr<R> randGenIn, double p_in = 0.0, double sigma_in=0.0){
            reg = ds_create_register(n, p_in, sigma_in);
            ds_set_state(reg, 0, 1, 0);
        }

        //! Enables an error model sim to begin a moment;
        //! Does nothing here
        void moment(){};

        //! Performs a controlled-X gate, (CNOT gate).
        //! \param control The control qubit id.
        //! \param target The target qubit id.
        void cx(int control, int target){
            //std::cout << "cx " << control << " "<< target << "\n";
            ds_cnot(reg, control, target, enable_errors);
        }

        //! Performs a controlled-Y gate -- NOT IMPLEMENTED
        //! \param control The control qubit id.
        //! \param target The target qubit id.
        //TODO: Implement
        //void cy(int control, int target);

        //! Performs a controlled-Z gate
        //! \param control The control qubit id.
        //! \param target The target qubit id.
        void cz(int control, int target){
            ds_Hadamard(reg, target, 0);
            ds_cnot(reg, control, target, enable_errors);
            ds_Hadamard(reg, target, 0);
        }

        //! Performs a CZ(theta) controlled rotation on control and target
        void cz(int control, int target, double theta){
            ds_cphase(reg, control, target, theta,  enable_errors);
        }

        //! Performs a Hadamard Gate
        //! \param qubit the qubit to perform the gate on.
        void h(int qubit){
            //std::cout << "h " << qubit << "\n";
            ds_Hadamard(reg, qubit, enable_errors);
        }

        //! Performs a Phase Gate
        //! \param qubit the qubit to perform the gate on.
        void s(int qubit){
            ds_zrot(reg, qubit, PI_2, enable_errors);
        }

        //! Performs an inverse Phase Gate
        //! \param qubit the qubit to perform the gate on.
        void sdg(int qubit){
            ds_zrot(reg, qubit, -PI_2, enable_errors);
        }

        //! Performs an X-Pauli gate
        //! \param qubit the qubit to perform the gate on.
        void x(int qubit){
            //std::cout << "x " << qubit << "\n";
            ds_X(reg, qubit, enable_errors);
        }

        //! Performs an Y-Pauli gate
        //! \param qubit the qubit to perform the gate on.
        void y(int qubit){
            ds_XZ(reg, qubit, enable_errors);
        }

        //! Performs an Z-Pauli gate
        //! \param qubit the qubit to perform the gate on.
        void z(int qubit){
            ds_Z(reg, qubit, enable_errors);
        }

        //! Performs an arbitrary X-rotation
        //! \param qubit the qubit to perform the gate on.
        //! \param theta the angle of rotation in radians.
        void x(int qubit, double theta){
            ds_xrot(reg, qubit, theta, enable_errors);
        }

        //! Performs an arbitrary Y-rotation
        //! \param qubit the qubit to perform the gate on.
        //! \param theta the angle of rotation in radians.
        void y(int qubit, double theta){
            ds_yrot(reg, qubit, theta, enable_errors);
        }

        //! Performs an arbitrary Z-rotation
        //! \param qubit the qubit to perform the gate on.
        //! \param theta the angle of rotation in radians.
        void z(int qubit, double theta){
            ds_zrot(reg, qubit, theta, enable_errors);
        }

        //! Performs an identity gate on a given qubit. In the base simulator this does nothing.
        //! \param qubit the qubit to perform the gate on.
        void idle(int qubit){}

        template<typename T>
        double getParityOneMeasurementProbability(T qubits){
            boost::container::small_vector<qubit_id_t, 20> temp{};
            for(auto q: qubits){
                temp.push_back(q);
            }
            return ds_set_measure_parity(reg, temp.size(), temp.data());
        }

        //! Measures a qubit in the Z basis.
        //! \param qubit
        //! \return The measurement outcome as a measurement_state;
        measurement_state measure(qubit_id_t qubit){
            measurement_state ret = ds_measure(reg, 1, &qubit) ? measurement_state::RANDOM_TRUE : measurement_state::RANDOM_FALSE;
            //std::cout << "m " << qubit << " "<< ((ret == measurement_state::RANDOM_TRUE)? 1:0) << "\n";
            return ret;
        }

        template<typename T>
        boost::container::small_vector<measurement_state, 20> measure_all(T qubits){
            boost::container::small_vector<qubit_id_t, 20> temp{};
            for(auto  q: qubits){
                temp.push_back(q);
            }
            auto measurement = ds_measure(reg, temp.size(), temp.data());
            boost::container::small_vector<measurement_state, 20> retval{};
            for (size_t i = 0; i < qubits.size(); ++i){
                retval.push_back((measurement & (1<<(qubits.size() - i - 1)))? measurement_state::RANDOM_TRUE : measurement_state::RANDOM_FALSE);
            }
            return retval;
        }


        template<typename T>
        boost::container::small_vector<measurement_state, 20> measure_reset_all(T qubits){
            auto retval = measure_all(qubits);

            for (size_t i = 0; i < retval.size(); i++) {
                if (retval[i] == measurement_state::RANDOM_TRUE || retval[i] == measurement_state::DETERMINISTIC_TRUE) {
                    x(qubits[i]);
                }
            }
            return retval;
        }

        template<typename T, typename U>
        double measure_all_set(T qubits, U required){
            int value = 0;
            for (size_t i = 0; i < qubits.size(); ++i){
                if((required[i] == measurement_state::DETERMINISTIC_TRUE) || (required[i] == measurement_state::RANDOM_TRUE)){
                    value |= (1<<(qubits.size() - i - 1));
                }
            }

            boost::container::small_vector<qubit_id_t, 20> temp{};
            for(auto q: qubits){
                temp.push_back(q);
            }

            return ds_set_measure(reg, temp.size(), temp.data(), value);
        }

        template<typename T, typename U>
        double measure_reset_all_set(T qubits, U required){
            auto retval = measure_all_set(qubits, required);

            for (size_t i = 0; i < required.size(); i++) {
                if (required[i] == measurement_state::RANDOM_TRUE || required[i] == measurement_state::DETERMINISTIC_TRUE) {
                    x(qubits[i]);
                }
            }
            return retval;
        }

        //! Measures a qubit in the Z basis and then resets it to single_pauli_state::ZERO.
        //! \param qubit
        //! \return The measurement outcome as a measurement_state;
        measurement_state measure_and_reset(int qubit){
            auto measurement = measure(qubit);
            if ((measurement == measurement_state::DETERMINISTIC_TRUE) || (measurement == measurement_state::RANDOM_TRUE)) {
                x(qubit);
            }
            return measurement;
        }

        //! Measures a qubit in the given pauli basis and then resets it to single_pauli_state::ZERO.
        //! \param qubit
        //! \param the basis to measure in
        //! \return The measurement outcome as a measurement_state;
        measurement_state measure_in_basis_and_reset(int qubit, single_pauli_basis b) {
            switch(b){
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

        //! Measures a qubit in the Z basis and then resets it to the state s;
        //! \param qubit
        //! \return The measurement outcome as a measurement_state;
        measurement_state measure_and_reset_to_state(int qubit, single_pauli_state state){
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

        //! Measures a qubit in the given pauli basis and then resets it to single_pauli_state::ZERO.
        //! \param qubit
        //! \param the basis to measure in
        //! \return The measurement outcome as a measurement_state;
        measurement_state measure_in_basis_and_reset_to_state(int qubit, single_pauli_basis basis, single_pauli_state state){
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

        //! Converts the state to a std::string that represents the tableu of the state.
        explicit operator std::string(){
            return ds_string(reg);
        };

        inline double getInnerProduct(const StateSimulator& other){
            return ds_inner_product(reg, other.reg);
        }

        //! compares equality the inner product
        inline bool operator==(const StateSimulator &other){
            return abs(1-ds_inner_product(reg, other.reg)) < (1.0E-200); // a very small number.
        }

        //! Compares if the stabiliziser parts are not equal.
        inline bool operator!=(const StateSimulator &other){
            return !((*this) == other);
        }

        //! Helper function that performs a rotation of a given qubit from the ZERO state to the given state;
        //! \param qubit - The qubit to rotate
        //! \param state - The state to rotate to
        void rotate_qubit(int qubit, single_pauli_state state){
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

        // void normalize(void)
        // void gaussian(void)
        // void ut_gaussian(void)
    };
}

#endif