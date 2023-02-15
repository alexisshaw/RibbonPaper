//
// Created by 13114347 on 25/05/2020.
//
#include <iostream>
#include <random>
#include "sim/chp.h"
#include <assert.h>

using namespace std;

bool test_pauli_product_phase();
bool test_identity();
bool test_bit_flip();
bool test_identity_2();
bool test_bit_flip_2();
bool test_epr();
bool test_phase_kickback_consume_s_state();
bool test_x();
bool test_phase_kickback_preserve_s_state();
bool test_kickback_vs_stabilizer();

std::shared_ptr<std::mt19937_64> getRandGen();
using namespace CodeTools;

int main(){
    cout << "Testing pauli_Product_phase(): "<< test_pauli_product_phase() << endl;
    cout << "Testing Identity: " << test_identity() << endl;
    cout << "Testing Bit Flip: " << test_bit_flip() << endl;
    cout << "Testing Identity on 2 qubits: " << test_identity_2() << endl;
    cout << "Testing Bit Flip with 2 qubits: " << test_bit_flip_2() << endl;
    cout << "Testing Einstein-Pauil-Rosen violation: " << test_epr() << endl;
    cout << "Testing Phase Kickback consuming an s state: " << test_phase_kickback_consume_s_state() << endl;
    cout << "Testing Bit Flip using an x gate: " << test_x() << endl;
    cout << "Testing Phase Kickback preserving the s state: " << test_phase_kickback_preserve_s_state() << endl;
    cout << "Testing for the stabilizers during kickback: " << test_kickback_vs_stabilizer() << endl;
    // TODO: Distillation Tests;

    return 0;
}

bool test_pauli_product_phase(){
    bool paulis[4][2] = {{false, false},
                         {true, false},
                         {true, true},
                         {false, true}};
    int expected[4][4] = {{0,0,0,0},{0,0,1,-1},{0,-1,0,1}, {0,1,-1,0}};

    assert(ChpSimulator<1>::pauli_product_phase(true, false, true, true) == 1);
    for (int i =0; i!= 4; i++) {
        for (int j = 0; j != 4; j++) {
            if (ChpSimulator<1>::pauli_product_phase(paulis[i][0], paulis[i][1], paulis[j][0], paulis[j][1]) != expected[i][j]) {
                return false;
            }
        }
    }
    return true;
}

bool test_identity(){
    auto s = ChpSimulator<1>(getRandGen());
    return s.measure_and_reset(0) == measurement_state::DETERMINISTIC_FALSE;
}

bool test_bit_flip(){
    auto s = ChpSimulator<1>(getRandGen());
    s.h(0);
    s.s(0);
    s.s(0);
    s.h(0);
    return s.measure_and_reset(0) == measurement_state::DETERMINISTIC_TRUE;
}
bool test_x(){
    auto s = ChpSimulator<1>(getRandGen());
    s.x(0);
    return s.measure_and_reset(0) == measurement_state::DETERMINISTIC_TRUE;
}

bool test_identity_2(){
    auto s = ChpSimulator<2>(getRandGen());
    return (s.measure_and_reset(0) == measurement_state::DETERMINISTIC_FALSE) &&
           (s.measure_and_reset(1) == measurement_state::DETERMINISTIC_FALSE);
}

bool test_bit_flip_2(){
    auto s = ChpSimulator<2>(getRandGen());
    s.h(0);
    s.s(0);
    s.s(0);
    s.h(0);
    return (s.measure_and_reset(0) == measurement_state::DETERMINISTIC_TRUE) &&
           (s.measure_and_reset(1) == measurement_state::DETERMINISTIC_FALSE);
}

bool test_epr(){
    auto s = ChpSimulator<2>(getRandGen());
    s.h(0);
    s.cx(0,1);
    auto v1 = (unsigned int)(s.measure(0));
    auto v2 = (unsigned int)(s.measure(1));
    return ((v1 & 2) != 2) && ((v2 & 2) == 2) && ((v1 & 1) == (v2 & 1));
}

bool test_phase_kickback_consume_s_state() {
    auto s = ChpSimulator<2>(getRandGen());
    s.h(1);
    s.s(1);
    s.h(0);
    s.cx(0,1);
    auto v1 = (unsigned int)(s.measure(1));
    if (v1 & 1){
        s.z(0);
    }
    s.s(0);
    s.h(0);
    auto v2 = (unsigned int)(s.measure(0));
    return ((v1 & 2) == 0) && (v2 == 3);
}

bool test_phase_kickback_preserve_s_state(){
    auto s = ChpSimulator<2>(getRandGen());

    // Prepare S state
    s.h(1);
    s.s(1);

    // Prepare test input.
    s.h(0);

    // Kickback
    s.cx(0,1);
    s.h(1);
    s.cx(0,1);
    s.h(1);

    // check
    s.s(0);
    s.h(0);
    auto m1 = unsigned(s.measure(0));
    s.s(1);
    s.h(1);
    auto m2 = unsigned(s.measure(1));
    return (m1 == 3) && (m2 == 3);
}

class FixedRandom{
public:
    int value;
    FixedRandom() = default;
    int operator()(){return value;};
    void seed(int seed){
    };
    static constexpr int min(){
        return 0;
    }
    static constexpr int max(){
        return 1;
    }
};

std::shared_ptr<FixedRandom> getFixedRandom(int value){
    std::shared_ptr<FixedRandom> gen = make_shared<FixedRandom>();
    gen->value = value;
    return gen;
}

bool test_kickback_vs_stabilizer(){
    auto s = ChpSimulator<3, FixedRandom>(getFixedRandom(0));

    s.h(2);
    s.cx(2,0);
    s.cx(2,1);
    s.s(0);
    s.s(1);
    s.h(0);
    s.h(1);
    s.h(2);
    const string c0 = "-Y..\n"
                      "-.Y.\n"
                      "+..X\n"
                      "----\n"
                      "+X.X\n"
                      "+.XX\n"
                      "+YYZ\n";
    bool comp0 = (string(s).compare(c0) == 0);
    auto v0 = unsigned(s.measure(0));
    const string c1 = "+X.X\n"
                      "-.Y.\n"
                      "+..X\n"
                      "----\n"
                      "+Z..\n"
                      "+.XX\n"
                      "+ZYY\n";

    bool comp1 = (string(s).compare(c1) == 0);
    auto v1 = unsigned(s.measure(1));
    const string c2 = "+X.X\n"
                      "+.XX\n"
                      "+..X\n"
                      "----\n"
                      "+Z..\n"
                      "+.Z.\n"
                      "-ZZZ\n";

    int comp2 = (string(s).compare(c2) == 0);
    auto v2 = unsigned(s.measure(2));
    int comp3 = (string(s) == c2);
    return (v0 == 0) && (v1 == 0) && (v2 == 3) && comp0 && comp1 && comp2 && comp3;
}

std::shared_ptr<std::mt19937_64> getRandGen(){
    std::random_device rd{"/dev/random"};
    std::seed_seq seq{rd(), rd(), rd(), rd()};
    std::shared_ptr<std::mt19937_64> gen = make_shared<std::mt19937_64>(seq);
    return gen;
}
