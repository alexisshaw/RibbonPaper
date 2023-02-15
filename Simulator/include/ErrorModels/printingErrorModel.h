//
// Created by 13114347 on 27/05/2020.
//

#ifndef STEANE_SIM_PRINTINGERRORMODEL_HPP
#define STEANE_SIM_PRINTINGERRORMODEL_HPP

#include <utility>
#include <iostream>

#include "magic_enum.h"
#include "BaseTypes.h"

namespace CodeTools {
    struct PrintingErrorModel {
        PrintingErrorModel() {
        }

        single_qubit_pauli h(int qubit) {
            std::cout << "h " << qubit << std::endl;
            return I;
        };

        single_qubit_pauli idle(int qubit) {
            std::cout << "idle " << qubit << std::endl;
            return I;
        };

        two_qubit_pauli cx(int q1, int q2) {
            std::cout << "cx " << q1 << " , " << q2 << std::endl;
            return std::make_pair(single_qubit_pauli::I, single_qubit_pauli::I);
        };

        two_qubit_pauli cz(int q1, int q2) {
            std::cout << "cz " << q1 << " , " << q2 << std::endl;
            return std::make_pair(single_qubit_pauli::I, single_qubit_pauli::I);
        };

        std::pair<single_qubit_pauli, single_qubit_pauli> measure_and_reset(int qubit) {
            std::cout << "m " << qubit << std::endl;
            return std::make_pair(single_qubit_pauli::I, single_qubit_pauli::I);
        };

        std::pair<single_qubit_pauli, single_qubit_pauli> measure_in_basis_and_reset(int qubit, single_pauli_basis b){
            std::cout << "m " << magic_enum::enum_name(b) << " " << qubit << std::endl;
            return std::make_pair(single_qubit_pauli::I, single_qubit_pauli::I);
        }

        std::pair<single_qubit_pauli, single_qubit_pauli> measure_in_basis_and_reset_to_state(int qubit, single_pauli_basis b, single_pauli_state state){
            std::cout << "m" << magic_enum::enum_name(b) << " " << qubit << " r " << agic_enum::enum_name(state) << std::endl;
            return std::make_pair(single_qubit_pauli::I, single_qubit_pauli::I);
        }
    };
};

#endif //STEANE_SIM_PRINTINGERRORMODEL_HPP
