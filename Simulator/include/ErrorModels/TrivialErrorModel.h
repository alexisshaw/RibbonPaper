//
// Created by 13114347 on 2/06/2020.
//

#ifndef STEANE_SIM_TRIVIALERRORMODEL_H
#define STEANE_SIM_TRIVIALERRORMODEL_H

#include <utility>
#include "BaseTypes.h"

namespace CodeTools{
    class TrivialErrorModel {
    public:
        static constexpr single_qubit_pauli h(int qubit) {
            return single_qubit_pauli::I;
        };

        static constexpr single_qubit_pauli idle(int qubit) {
            return single_qubit_pauli::I;
        };

        static constexpr two_qubit_pauli cx(int q1, int q2) {
            return std::make_pair(single_qubit_pauli::I, single_qubit_pauli::I);
        };

        static constexpr two_qubit_pauli cz(int q1, int q2) {
            return std::make_pair(single_qubit_pauli::I, single_qubit_pauli::I);
        };

        static constexpr std::pair<single_qubit_pauli, single_qubit_pauli> measure_and_reset(int qubit) {
            return std::make_pair(single_qubit_pauli::I, single_qubit_pauli::I);
        };

        static constexpr std::pair<single_qubit_pauli, single_qubit_pauli> measure_in_basis_and_reset(int qubit, single_pauli_basis b){
            return std::make_pair(single_qubit_pauli::I, single_qubit_pauli::I);
        }

        static constexpr std::pair<single_qubit_pauli, single_qubit_pauli> measure_in_basis_and_reset_to_state(int qubit, single_pauli_basis b, single_pauli_state state){
            return std::make_pair(single_qubit_pauli::I, single_qubit_pauli::I);
        }
    };
}
#endif //STEANE_SIM_TRIVIALERRORMODEL_H
