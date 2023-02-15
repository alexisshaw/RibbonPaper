//
// Created by 13114347 on 26/05/2020.
//

#ifndef STEANE_SIM_BALANCEDERRORMODEL_H
#define STEANE_SIM_BALANCEDERRORMODEL_H

#include <random>
#include "codes/Steane.h"

namespace CodeTools {
    struct BalancedErrorModel {
        using R = std::mt19937_64;
        R rand_gen;

        typename R::result_type const p_cx = 0;
        typename R::result_type const p_cz = 0;
        typename R::result_type const p_h = 0;
        typename R::result_type const p_idle = 0;
        typename R::result_type const p_measure = 0;


        BalancedErrorModel() {
            std::random_device rd{"/dev/random"};
            std::seed_seq seq{rd(), rd(), rd(), rd()};
            rand_gen = R(seq);
        }

        BalancedErrorModel(long double p_cx_in,long double p_h_in,long double p_idle_in, long double p_measure_in):
            p_cx{static_cast<R::result_type>((((p_cx_in * 16.0L)/15.0L) * static_cast<long double>((R::max() - R::min())))) + R::min()},
            p_cz{static_cast<R::result_type>((((p_cx_in * 16.0L)/15.0L) * static_cast<long double>((R::max() - R::min())))) + R::min()},
            p_h{static_cast<R::result_type>((((p_h_in * 4.0L)/3.0L) * static_cast<long double>((R::max() - R::min())))) + R::min()},
            p_idle{static_cast<R::result_type>((((p_idle_in * 4.0L)/3.0L) * static_cast<long double>((R::max() - R::min())))) + R::min()},
            p_measure{static_cast<R::result_type>((((p_measure_in * 4.0L)/3.0L) * static_cast<long double>((R::max() - R::min())))) + R::min()}{
            std::random_device rd("/dev/random");
            std::seed_seq seq{rd(), rd(), rd(), rd()};
            rand_gen = R{seq};
        }

        single_qubit_pauli h(int qubit){
            if(rand_gen() < p_h){
                R::result_type rn = rand_gen();
                return single_qubit_pauli(rn & 3);
            }
            return single_qubit_pauli::I;
        };

        single_qubit_pauli idle(int qubit){
            if(rand_gen() < p_idle){
                R::result_type rn = rand_gen();
                return single_qubit_pauli(rn & 3);
            }
            return single_qubit_pauli::I;
        };

        two_qubit_pauli cx(int q1, int q2){
            if(rand_gen() < p_cx){
                R::result_type rn = rand_gen();
                return std::make_pair(single_qubit_pauli(rn & 3), single_qubit_pauli((rn >> 2) & 3));
            }
            return std::make_pair(single_qubit_pauli::I, single_qubit_pauli::I);
        };

        two_qubit_pauli cz(int q1, int q2){
            if(rand_gen() < p_cz){
                R::result_type rn = rand_gen();
                return std::make_pair(single_qubit_pauli(rn & 3), single_qubit_pauli((rn >> 2) & 3));
            }
            return std::make_pair(single_qubit_pauli::I, single_qubit_pauli::I);
        };

        std::pair<single_qubit_pauli, single_qubit_pauli> measure_and_reset(int qubit){
            if(rand_gen() < p_measure){
                R::result_type rn = rand_gen();
                return std::make_pair(single_qubit_pauli(rn & 3), single_qubit_pauli::I);
            }
            return std::make_pair(single_qubit_pauli::I, single_qubit_pauli::I);
        };

        std::pair<single_qubit_pauli, single_qubit_pauli> measure_in_basis_and_reset(int qubit, single_pauli_basis b){
            if(rand_gen() < p_measure){
                R::result_type rn = rand_gen();
                return std::make_pair(single_qubit_pauli(rn & 3), single_qubit_pauli::I);
            }
            return std::make_pair(single_qubit_pauli::I, single_qubit_pauli::I);
        }

        std::pair<single_qubit_pauli, single_qubit_pauli> measure_in_basis_and_reset_to_state(int qubit, single_pauli_basis b, single_pauli_state state){
            if(rand_gen() < p_measure){
                R::result_type rn = rand_gen();
                return std::make_pair(single_qubit_pauli(rn & 3), single_qubit_pauli::I);
            }
            return std::make_pair(single_qubit_pauli::I, single_qubit_pauli::I);
        }

    };
}

#endif //STEANE_SIM_BALANCEDERRORMODEL_H