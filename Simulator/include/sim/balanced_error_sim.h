//
// Created by 13114347 on 24/06/2020.
//

#ifndef STEANE_SIM_BALANCED_ERROR_SIM_H
#define STEANE_SIM_BALANCED_ERROR_SIM_H

#include <string>
#include <random>
#include <memory>
#include "BaseTypes.h"

namespace CodeTools{
//! An adaptor class that adds balanced Random Errors to an underlying sim class
template<typename BaseSim, typename R=std::mt19937_64>
class BalancedErrorSim {
public:
    BaseSim baseSim;
    static constexpr int N{BaseSim::N};
private:
    std::shared_ptr<R> rand_gen;

    using RandomValue = typename R::result_type;
    RandomValue const p_single = 0;
    RandomValue const p_dual = 0;
    RandomValue const p_measure_pre = 0;
    RandomValue const p_measure_post = 0;

    std::bitset<N> tracker{};

    void maybe_single_error(int qubit){
        if((*rand_gen)() < p_single) {
            RandomValue rn = (*rand_gen)();
            baseSim.inject_single_pauli_error(qubit, single_qubit_pauli(rn & 3));
        }
        tracker.set(qubit);
    }
    void maybe_dual_error(int q1, int q2){
        if((*rand_gen)() < p_dual){
            RandomValue rn = (*rand_gen)();
            baseSim.inject_single_pauli_error(q1, single_qubit_pauli(rn & 3));
            baseSim.inject_single_pauli_error(q2, single_qubit_pauli((rn >> 2) & 3));
        }
        tracker.set(q1);
        tracker.set(q2);
    }

    void maybe_measure_error_pre(int qubit){
        if((*rand_gen)() < p_measure_pre) {
            RandomValue rn = (*rand_gen)();
            baseSim.inject_single_pauli_error(qubit, single_qubit_pauli(rn & 3));
        }
    }
    void maybe_measure_error_post(int qubit){
        if((*rand_gen)() < p_measure_post) {
            RandomValue rn = (*rand_gen)();
            baseSim.inject_single_pauli_error(qubit, single_qubit_pauli(rn & 3));
        }
        tracker.set(qubit);
    }
public:
    //! A base simulator class, operations on which occur error free
    using error_free_sim = typename BaseSim::error_free_sim;

    BalancedErrorSim(const BaseSim& b, std::shared_ptr<R> rand_gen_in): baseSim{b}, rand_gen(rand_gen_in){
    }

    BalancedErrorSim(long double p_single_in,
                     long double p_dual_in,
                     long double p_measure_pre_in,
                     long double p_measure_post_in,
                     const BaseSim& b,
                     std::shared_ptr<R> rand_gen_in):
            baseSim{b}, rand_gen{rand_gen_in},
                     p_single(static_cast<RandomValue>((p_single_in * 4.0L)/3.0L*static_cast<long double>((R::max() - R::min()))) + R::min()),
    p_dual(static_cast<RandomValue>((p_dual_in * 16.0L)/15.0L*static_cast<long double>((R::max() - R::min()))) + R::min()),
    p_measure_pre(static_cast<RandomValue>((p_measure_pre_in * 4.0L)/3.0L*static_cast<long double>((R::max() - R::min()))) + R::min()),
    p_measure_post(static_cast<RandomValue>((p_measure_post_in * 4.0L)/3.0L*static_cast<long double>((R::max() - R::min()))) + R::min())
{}

    void inject_single_pauli_error(int qubit, single_qubit_pauli p){
        baseSim.inject_single_pauli_error(qubit,p);
    };

    //! Enables an error model sim to begin a moment;
    //! Does nothing here
    void moment(){
        for(int i=0; i < N; ++i){
            if(!tracker.test(i)){
                idle(i);
            }
        }
        tracker.reset();

    };

    //! Performs a controlled-X gate, (CNOT gate).
    //! \param control The control qubit id.
    //! \param target The target qubit id.
    void cx(int control, int target){
        maybe_dual_error(control, target);
        baseSim.cx(control, target);
    }

    //! Performs a controlled-Y gate
    //! \param control The control qubit id.
    //! \param target The target qubit id.
    void cy(int control, int target){
        maybe_dual_error(control, target);
        baseSim.cy(control,target);
    }

    //! Performs a controlled-Z gate
    //! \param control The control qubit id.
    //! \param target The target qubit id.
    void cz(int control, int target){
        maybe_dual_error(control, target);
        baseSim.cz(control, target);
    }

    //! Performs a Hadamard Gate
    //! \param qubit the qubit to perform the gate on.
    void h(int qubit){
        maybe_single_error(qubit);
        baseSim.h(qubit);
    }

    //! Performs a Phase Gate
    //! \param qubit the qubit to perform the gate on.
    void s(int qubit){
        maybe_single_error(qubit);
        baseSim.s(qubit);
    }

    //! Performs an inverse Phase Gate
    //! \param qubit the qubit to perform the gate on.
    void sdg(int qubit){
        maybe_single_error(qubit);
        baseSim.sdg(qubit);
    }

    //! Performs an X-Pauli gate
    //! \param qubit the qubit to perform the gate on.
    void x(int qubit){
        maybe_single_error(qubit);
        baseSim.x(qubit);
    }

    //! Performs an Y-Pauli gate
    //! \param qubit the qubit to perform the gate on.
    void y(int qubit){
        maybe_single_error(qubit);
        baseSim.y(qubit);
    }

    //! Performs an Z-Pauli gate
    //! \param qubit the qubit to perform the gate on.
    void z(int qubit){
        maybe_single_error(qubit);
        baseSim.z(qubit);
    }

    //! Performs an identity gate on a given qubit.
    //! \param qubit the qubit to perform the gate on.
    void idle(int qubit){
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
    boost::container::small_vector<measurement_state, 20> measure_all(T qubits){
        boost::container::small_vector<measurement_state, 20> retval{};
        for (auto qubit: qubits){
             maybe_measure_error_pre(qubit);
        }
        for (auto qubit: qubits){
            retval.push_back(measure(qubit));
        }
        for (auto qubit: qubits){
            maybe_measure_error_post(qubit);
        }
        return retval;
    }


    //! Measures a qubit in the Z basis and then resets it to single_pauli_state::ZERO.
    //! \param qubit
    //! \return The measurement outcome as a measurement_state;
    measurement_state measure_and_reset(int qubit){
        maybe_measure_error_pre(qubit);
        auto result = baseSim.measure_and_reset(qubit);
        maybe_measure_error_post(qubit);
        return result;
    }

    //! Measures a qubit in the given pauli basis and then resets it to single_pauli_state::ZERO.
    //! \param qubit
    //! \param the basis to measure in
    //! \return The measurement outcome as a measurement_state;
    measurement_state measure_in_basis_and_reset(int qubit, single_pauli_basis b){
        maybe_measure_error_pre(qubit);
        auto result = baseSim.measure_in_basis_and_reset(qubit, b);
        maybe_measure_error_post(qubit);
        return result;
    }

    //! Measures a qubit in the Z basis and then resets it to the state s;
    //! \param qubit
    //! \return The measurement outcome as a measurement_state;
    measurement_state measure_and_reset_to_state(int qubit, single_pauli_state s){
        maybe_measure_error_pre(qubit);
        auto result = baseSim.measure_and_reset_to_state(qubit, s);
        maybe_measure_error_post(qubit);
        return result;
    }

    //! Measures a qubit in the given pauli basis and then resets it to single_pauli_state::ZERO.
    //! \param qubit
    //! \param the basis to measure in
    //! \return The measurement outcome as a measurement_state;
    measurement_state measure_in_basis_and_reset_to_state(int qubit, single_pauli_basis b, single_pauli_state s){
        maybe_measure_error_pre(qubit);
        auto result = baseSim.measure_in_basis_and_reset_to_state(qubit, b, s);
        maybe_measure_error_post(qubit);
        return result;
    }

    //! Converts the state to a std::string that represents the tableu of the state.
    explicit operator std::string(){
        return std::string(baseSim);
    }

    //! compares equality using only the stabilizer part of the tableu, does not normalize in the comparison.
    //! in future we should change this to compare based on whether the tables represent the same states by gaussian
    //! elimination.
    inline bool operator==(BalancedErrorSim<BaseSim, R> &other){
        return baseSim == other.baseSim;
    }

    //! Compares if the stabiliziser parts are not equal.
    inline bool operator!=(BalancedErrorSim<BaseSim, R> &other){
        return !((*this) == other);
    }

    //! Helper function that performs a rotation of a given qubit from the ZERO state to the given state;
    //! This fuction does not add errors.
    //! \param qubit - The qubit to rotate
    //! \param state - The state to rotate to
    inline void rotate_qubit(int qubit, single_pauli_state state){
        return baseSim.rotate_qubit(qubit, state);
    }

    //! Function that returns an error free copy of the base class, to enable clean compsition of Error correction simulations.
    error_free_sim& get_error_free(){
        return baseSim.get_error_free();
    }

    //! Function that returns an error free copy of the base class, to enable clean compsition of Error correction simulations.
    const error_free_sim& get_error_free() const{
        return baseSim.get_error_free();
    }
};
}

#endif //STEANE_SIM_BALANCED_ERROR_SIM_H
