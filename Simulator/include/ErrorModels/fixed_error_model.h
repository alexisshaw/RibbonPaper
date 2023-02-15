//
// Created by 13114347 on 13/06/2020.
//

#include <utility>
#include "BaseTypes.h"

#ifndef STEANE_SIM_FIXED_ERROR_MODEL_HPP
#define STEANE_SIM_FIXED_ERROR_MODEL_HPP

namespace CodeTools{
    class FixedErrorModel {
    public:
        single_qubit_pauli single_error = single_qubit_pauli::I;
        two_qubit_pauli dual_error = std::make_pair(single_qubit_pauli::I, single_qubit_pauli::I);
        int single_count = 0;
        int dual_count = 0;
        int single_count_target = -1;
        int dual_count_target = -1;

    public:
        FixedErrorModel(single_qubit_pauli p, int target): single_error(p), single_count_target(target){};
        FixedErrorModel(two_qubit_pauli p, int target): dual_error(std::move(p)), dual_count_target(target){};
        FixedErrorModel() = default;

        single_qubit_pauli h(int qubit) {
            if(single_count_target == (single_count++)){
                //std::cout << "h " << qubit << std::endl;
                return single_error;
            } else {
                return single_qubit_pauli::I;
            }
        }

        single_qubit_pauli idle(int qubit) {
            if(single_count_target == (single_count++)){
                //std::cout << "idle " << qubit << std::endl;
                return single_error;
            } else {
                return single_qubit_pauli::I;
            }
        }

        two_qubit_pauli cx(int q1, int q2) {
            if(dual_count_target == (dual_count++)){
                //std::cout << "cx " << q1 << " , " << q2 << std::endl;
                return dual_error;
            } else {
                return std::make_pair(single_qubit_pauli::I, single_qubit_pauli::I);
            }
        }

        two_qubit_pauli cz(int q1, int q2) {
            if(dual_count_target == (dual_count++)){
                //std::cout << "cz " << q1 << " , " << q2 << std::endl;
                return dual_error;
            } else {
                return std::make_pair(single_qubit_pauli::I, single_qubit_pauli::I);
            }
        }

        std::pair<single_qubit_pauli, single_qubit_pauli> measure_and_reset(int qubit) {
            single_qubit_pauli first;
            single_qubit_pauli second;
            if(single_count_target == (single_count++)){
                //std::cout << "mr_pre " << qubit << std::endl;
                first = single_error;
            } else {
                first = single_qubit_pauli::I;
            }
            if(single_count_target == (single_count++)){
                //std::cout << "mr_post " << qubit << std::endl;
                second = single_error;
            } else {
                second = single_qubit_pauli::I;
            }
            return std::make_pair(first, second);
        };

        std::pair<single_qubit_pauli, single_qubit_pauli> measure_in_basis_and_reset(int qubit, single_pauli_basis b){
            single_qubit_pauli first;
            single_qubit_pauli second;
            if(single_count_target == (single_count++)){
                //std::cout << "mbr_pre " << qubit << std::endl;
                first = single_error;
            } else {
                first = single_qubit_pauli::I;
            }
            if(single_count_target == (single_count++)){
                //std::cout << "mbr_post " << qubit << std::endl;
                second = single_error;
            } else {
                second = single_qubit_pauli::I;
            }
            return std::make_pair(first, second);
        }

        std::pair<single_qubit_pauli, single_qubit_pauli> measure_in_basis_and_reset_to_state(int qubit, single_pauli_basis b, single_pauli_state state){
            single_qubit_pauli first;
            single_qubit_pauli second;
            if(single_count_target == (single_count++)){
                //std::cout << "mbs_pre " << qubit << std::endl;
                first = single_error;
            } else {
                first = single_qubit_pauli::I;
            }
            if(single_count_target == (single_count++)){
                //std::cout << "mbs_post " << qubit << std::endl;
                second = single_error;
            } else {
                second = single_qubit_pauli::I;
            }
            return std::make_pair(first, second);
        }
    };
}

#endif //STEANE_SIM_FIXED_ERROR_MODEL_HPP
