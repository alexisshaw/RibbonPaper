//
// Created by 13114347 on 8/10/2020.
//

#ifndef STEANE_SIM_UTILITYFUNCTIONS_H
#define STEANE_SIM_UTILITYFUNCTIONS_H
#include <bitset>
#include <iostream>
#include <variant>
#include "BaseTypes.h"

namespace CorrectorFactory {
    using CodeTools::error_correction_list;
    using CodeTools::single_qubit_pauli;

    template<int qubits>
    std::bitset<2*qubits> to_bitset(const error_correction_list& l){
        std::bitset<2*qubits> b;
        for (const auto& pauli: l){
            switch(pauli.second){
                case single_qubit_pauli::X:
                    b.set(pauli.first + qubits);
                    break;
                case single_qubit_pauli::Z:
                    b.set(pauli.first);
                    break;
                case single_qubit_pauli::Y:
                    b.set(pauli.first + qubits);
                    b.set(pauli.first);
                    break;
                default:
                    break;
            }
        }
        return b;
    }

    template<int qubits>
    error_correction_list to_correction_list(std::bitset<2*qubits> b){
        error_correction_list l{};
        for(size_t i=0; i< qubits; ++i) {
            if (b.test(i) && b.test(i + qubits)) {
                l.push_back(std::make_pair(i, single_qubit_pauli::Y));
            } else if (b.test(i)) {
                l.push_back(std::make_pair(i, single_qubit_pauli::Z));
            } else if (b.test(i+qubits)){
                l.push_back(std::make_pair(i, single_qubit_pauli::X));
            }

        }
        return l;
    }

    template<typename T>
    bool compare_ErrorSyndrome(const error_correction_list& A, const error_correction_list& B){
        if(A == B){
            return true;
        }
        constexpr size_t qubits = T::REQUIRED_QUBITS + T::REQUIRED_ANCS;
        std::bitset<2*qubits> A_bitset(to_bitset<qubits>(A));
        std::bitset<2*qubits> B_bitset(to_bitset<qubits>(B));
        A_bitset = A_bitset ^ B_bitset;

        for(const auto &stab: T::STABS){
            if(A_bitset.test(stab.front()))
                for (int qubit_no: stab)
                    A_bitset.flip(qubit_no);
            if(A_bitset.test(stab.front() + qubits))
                for (int qubit_no : stab)
                    A_bitset.flip(qubit_no + qubits);
        }
        return A_bitset.none();
    }

    // This function tries to find a lower weight correction that is equivalent to the given correction list,
    // In this version it cheats by only comparing with the stabilisers, in order to try and get faster.
    // TODO: Make this run over the power set of the stab_bitsets, or do something more smart.
    template<typename T>
    error_correction_list simplifyErrorSyndrome(const error_correction_list& A){
        constexpr size_t qubits = T::REQUIRED_QUBITS + T::REQUIRED_ANCS;
        //Generate an array of bitsets.
        std::array<std::bitset<2*qubits>, 2*T::NUM_STABS> stab_bitsets{};
        for(int i=0; i< T::NUM_STABS; ++i){
            for(auto qubit_no: T::STABS[i]) {
                stab_bitsets[i].set(qubit_no);
                stab_bitsets[i+T::NUM_STABS].set(qubit_no + qubits);
            }
        }
        if(A.size() <= 1){
            return A;
        }

        std::bitset<2*qubits> A_bitset(to_bitset<qubits>(A));
        std::bitset<2*qubits> current_best(A_bitset);
        for(size_t i = 0; i < stab_bitsets.size(); i++){
            std::bitset<2*qubits> curr = A_bitset ^ stab_bitsets[i];
            if(curr.count() < current_best.count())
                current_best = curr;

        }

        std::cout << std::flush;
        return to_correction_list<qubits>(current_best);
    }



    template <typename Key, typename T>
    std::multimap<Key, T> simplify_multimap(const std::multimap<Key, T>& m){
        std::set<std::pair<const Key, T>> temp_set{m.begin(), m.end()};
        return std::multimap<Key, T>(temp_set.begin(), temp_set.end());
    }

}

//using two_error_key = std::pair<CodeTools::two_qubit_pauli, int>;
//using single_error_key = std::pair<CodeTools::single_qubit_pauli, int>;
//using error_key = std::variant<single_error_key, two_error_key>;

namespace std {
    template <typename TChar, typename TCharTraits>
    inline basic_ostream<TChar, TCharTraits>&
    operator<<(basic_ostream<TChar, TCharTraits> & stream, const CodeTools::single_qubit_pauli &p) {
        stream << magic_enum::enum_name(p);
        return stream;
    }

    template <typename TChar, typename TCharTraits>
    inline basic_ostream<TChar, TCharTraits>&
    operator<<(basic_ostream<TChar, TCharTraits> & stream, const CodeTools::single_pauli_state &s) {
        stream << magic_enum::enum_name(s);
        return stream;
    }

//    template <typename TChar, typename TCharTraits>
//    inline basic_ostream<TChar, TCharTraits>&
//    operator<<(basic_ostream<TChar, TCharTraits> & stream, const error_key &e){
//        std::visit([&stream](auto&& arg){
//            stream << arg;
//        }, e);
//        return stream;
//    }

    inline std::ostream &operator<<(std::ostream &out, const CodeTools::single_qubit_pauli &p) {
        out << magic_enum::enum_name(p);
        return out;
    }
}
#endif //STEANE_SIM_UTILITYFUNCTIONS_H
