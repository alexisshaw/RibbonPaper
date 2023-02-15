//
// Created by 13114347 on 24/06/2020.
//

#ifndef STEANE_SIM_AUTO_IDLER_H
#define STEANE_SIM_AUTO_IDLER_H

#include <string>
#include <random>
#include <memory>
#include "BaseTypes.h"

namespace CodeTools{
//! An adaptor class that adds balanced Random Errors to an underlying sim class
template<typename BaseSim>
class AutoIdlerSim {
public:
    BaseSim baseSim;
    static constexpr auto N{BaseSim::N};
private:
    std::array<int, N> steps{};

    void update(int q1, int q2){
        int &sq1 = steps[q1], &sq2 = steps[q2];
        int min_index = (sq1 < sq2)? q1: q2,  max_index = (sq1 < sq2)? q2: q1;

        while (steps[min_index] < steps[max_index]) {
            idle(min_index);
        }
    }

    void global_update(){
        int max_steps = *(std::max_element(begin(steps), end(steps)));

        for (size_t q=0; q<steps.size(); q++)
            while (steps[q] < max_steps)
                idle(static_cast<int>(q));
    }

public:
    //! A base simulator class, operations on which occur error free
    using error_free_sim = typename BaseSim::error_free_sim;

    AutoIdlerSim(const BaseSim& b):
            baseSim{b}{}

    AutoIdlerSim(BaseSim&& b):
            baseSim{std::move(b)}{}

    void inject_single_pauli_error(int qubit, single_qubit_pauli p){
        baseSim.inject_single_pauli_error(qubit,p);
    };

    //! Enables an error model sim to begin a moment;
    //! Does nothing here
    void moment(){
        global_update();
    };

    //! Performs a controlled-X gate, (CNOT gate).
    //! \param control The control qubit id.
    //! \param target The target qubit id.
    void cx(int control, int target){
        update(control, target);
        baseSim.cx(control, target);
        steps[control]++;
        steps[target]++;
    }

    //! Performs a controlled-Y gate
    //! \param control The control qubit id.
    //! \param target The target qubit id.
    void cy(int control, int target){
        update(control, target);
        baseSim.cy(control,target);
        steps[control]++;
        steps[target]++;
    }

    //! Performs a controlled-Z gate
    //! \param control The control qubit id.
    //! \param target The target qubit id.
    void cz(int control, int target){
        update(control, target);
        baseSim.cz(control, target);
        steps[control]++;
        steps[target]++;
    }

    //! Performs a controlled-Z gate
    //! \param control The control qubit id.
    //! \param target The target qubit id.
    void cz(int control, int target, double theta){
        update(control, target);
        baseSim.cz(control, target, theta);
        steps[control]++;
        steps[target]++;
    }

    //! Performs a Hadamard Gate
    //! \param qubit the qubit to perform the gate on.
    void h(int qubit){
        baseSim.h(qubit);
        steps[qubit]++;
    }

    //! Performs a Phase Gate
    //! \param qubit the qubit to perform the gate on.
    void s(int qubit){
        baseSim.s(qubit);
        steps[qubit]++;
    }

    //! Performs an inverse Phase Gate
    //! \param qubit the qubit to perform the gate on.
    void sdg(int qubit){
        baseSim.sdg(qubit);
        steps[qubit]++;
    }

    //! Performs an X-Pauli gate
    //! \param qubit the qubit to perform the gate on.
    void x(int qubit){
        baseSim.x(qubit);
        steps[qubit]++;
    }

    //! Performs an Y-Pauli gate
    //! \param qubit the qubit to perform the gate on.
    void y(int qubit){
        baseSim.y(qubit);
        steps[qubit]++;
    }

    //! Performs an Z-Pauli gate
    //! \param qubit the qubit to perform the gate on.
    void z(int qubit){
        baseSim.z(qubit);
        steps[qubit]++;
    }

    //! Performs an X-Pauli gate
    //! \param qubit the qubit to perform the gate on.
    void x(int qubit, double theta){
        baseSim.x(qubit, theta);
        steps[qubit]++;
    }

    //! Performs an Y-Pauli gate
    //! \param qubit the qubit to perform the gate on.
    void y(int qubit, double theta){
        baseSim.y(qubit, theta);
        steps[qubit]++;
    }

    //! Performs an Z-Pauli gate
    //! \param qubit the qubit to perform the gate on.
    void z(int qubit, double theta){
        baseSim.z(qubit,theta);
        steps[qubit]++;
    }

    //! Performs an identity gate on a given qubit.
    //! \param qubit the qubit to perform the gate on.
    void idle(int qubit){
        baseSim.idle(qubit);
        steps[qubit]++;
    }

    //! Measures a qubit in the Z basis
    //! \param qubit
    //! \return The measurement outcome as a measurement_state;
    measurement_state measure(int qubit){
        global_update();
        auto result = baseSim.measure(qubit);
        return result;
    }

    auto measure_all(auto qubits){
        boost::container::small_vector<measurement_state, 20> retval{};
        global_update();
        return baseSim.measure_all(qubits);
    }

    auto measure_reset_all(auto qubits){
        boost::container::small_vector<measurement_state, 20> retval{};
        global_update();
        return baseSim.measure_reset_all(qubits);
    }

    auto measure_reset_all_set(auto qubits, auto required){
        boost::container::small_vector<measurement_state, 20> retval{};
        // global_update();
        // for (auto qubit: qubits){
        //     retval.push_back(measure(qubit));
        // }
        return baseSim.measure_reset_all_set(qubits, required);
    }


    //! Measures a qubit in the Z basis and then resets it to single_pauli_state::ZERO.
    //! \param qubit
    //! \return The measurement outcome as a measurement_state;
    measurement_state measure_and_reset(int qubit){
        global_update();
        auto result = baseSim.measure_and_reset(qubit);
        return result;
    }

    //! Measures a qubit in the given pauli basis and then resets it to single_pauli_state::ZERO.
    //! \param qubit
    //! \param the basis to measure in
    //! \return The measurement outcome as a measurement_state;
    measurement_state measure_in_basis_and_reset(int qubit, single_pauli_basis b){
        global_update();
        auto result = baseSim.measure_in_basis_and_reset(qubit, b);
        return result;
    }

    //! Measures a qubit in the Z basis and then resets it to the state s;
    //! \param qubit
    //! \return The measurement outcome as a measurement_state;
    measurement_state measure_and_reset_to_state(int qubit, single_pauli_state s){
        global_update();
        auto result = baseSim.measure_and_reset_to_state(qubit, s);
        return result;
    }

    //! Measures a qubit in the given pauli basis and then resets it to single_pauli_state::ZERO.
    //! \param qubit
    //! \param the basis to measure in
    //! \return The measurement outcome as a measurement_state;
    measurement_state measure_in_basis_and_reset_to_state(int qubit, single_pauli_basis b, single_pauli_state s){
        global_update();
        auto result = baseSim.measure_in_basis_and_reset_to_state(qubit, b, s);
        return result;
    }

    //! Converts the state to a std::string that represents the tableu of the state.
    explicit operator std::string(){
        return std::string(baseSim);
    }

    //! compares equality using only the stabilizer part of the tableu, does not normalize in the comparison.
    //! in future we should change this to compare based on whether the tables represent the same states by gaussian
    //! elimination.
    inline bool operator==(AutoIdlerSim<BaseSim> &other){
        return baseSim == other.baseSim;
    }

    //! gets the inner product of two states if defined.
    inline double getInnerProduct(AutoIdlerSim<BaseSim> &other){
        return baseSim.getInnerProduct(other.baseSim);
    }

    inline auto getParityOneMeasurementProbability(auto qubits){
        return baseSim.getParityOneMeasurementProbability(qubits);
    }

    //! Compares if the stabiliziser parts are not equal.
    template<typename R>
    inline bool operator!=(AutoIdlerSim<BaseSim> &other){
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
