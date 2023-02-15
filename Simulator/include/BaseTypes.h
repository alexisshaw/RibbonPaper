//
// Created by 13114347 on 2/06/2020.
//

#ifndef STEANE_SIM_BASETYPES_H
#define STEANE_SIM_BASETYPES_H

#include <utility>
#include <vector>
#include <boost/container/static_vector.hpp>
#include <boost/container/small_vector.hpp>

namespace CodeTools {
    enum class single_pauli_state : uint8_t {
        ZERO, ONE, PLUS, MINUS, I_PLUS, I_MINUS
    };

    enum class single_qubit_pauli : uint8_t {
        I, Z, X, Y
    };

    enum class single_pauli_basis : uint8_t {
        Z, X, Y
    };

    enum class measurement_state : uint8_t {
        RANDOM_FALSE,
        RANDOM_TRUE,
        DETERMINISTIC_FALSE,
        DETERMINISTIC_TRUE
    };
    inline bool measurement_state_value(measurement_state in) {
        switch(in){
            case measurement_state::RANDOM_FALSE:
                return false;
            case measurement_state::RANDOM_TRUE:
                return true;
            case measurement_state::DETERMINISTIC_FALSE:
                return false;
            case measurement_state::DETERMINISTIC_TRUE:
                return true;
            default:
                return true;
        }
    }
    inline bool measurement_state_random(measurement_state in) {
        switch(in){
            case measurement_state::RANDOM_FALSE:
                return true;
            case measurement_state::RANDOM_TRUE:
                return true;
            case measurement_state::DETERMINISTIC_FALSE:
                return false;
            case measurement_state::DETERMINISTIC_TRUE:
                return false;
            default:
                return true;
        }
    }


    using two_qubit_pauli = std::pair<single_qubit_pauli, single_qubit_pauli>;

    using error_correction_list = boost::container::small_vector<std::pair<uint8_t, single_qubit_pauli>, 8>;

    template<size_t N>
    using static_error_correction_list = boost::container::static_vector<std::pair<uint8_t, single_qubit_pauli>, N>;
}

#endif //STEANE_SIM_BASETYPES_H
