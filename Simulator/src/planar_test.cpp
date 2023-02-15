#include <iostream>
#include "codes/Planar.h"
#include "sim/chp.h"
#include "sim/fixed_error_sim_post.h"
#include "sim/auto_idler.h"

#include "prettyprint.h"
#include "magic_enum.h"

// #include "util/FTCorrectorFactory.h"
// #include "util/testRunner.h"
// #include "absl/container/flat_hash_map.h"

#define DISTANCE 3

using namespace magic_enum::ostream_operators;
// using namespace boost::program_options;


//#include "absl/container/flat_hash_map.h"

// Takes a specific numered location in the lattice and distance d and specifies the relavant stabilizer indices, type = 0 for vertices and type = 1 for plaquettes

using namespace CodeTools;
using CodeTools::ChpSimulator;
using CodeTools::FixedErrorSimPost;
using CodeTools::AutoIdlerSim;
// using CodeTools::ErrorTracker;
// using CodeTools::TeeSim;
// using CodeTools::BalancedErrorSim;
// using CodeTools::single_pauli_state;


constexpr int distance = DISTANCE;
int limit = 100;
std::vector<double> p_series{};
constexpr auto initial_state = single_pauli_state::ZERO;
int nc = 40;
using PlanarBase = Planar::PlanarBase<distance>;
// using Sim = StateSimulator<PlanarBase::TOTAL_QUBITS>;
// using Planar_T = Planar::Planar<distance, Sim>;

template <class sim>
using planar_d = Planar::Planar<distance, sim>;
static constexpr int TOTAL_QUBITS  = PlanarBase::TOTAL_QUBITS;

auto get_rand_gen(){
    std::random_device rd{"/dev/urandom"};
    std::seed_seq seq{rd(), rd(), rd(), rd()};
    return std::make_shared<std::mt19937_64>(seq);
}

auto getPlanarErrorSim(auto rand_gen, auto errorOp, int error_n, single_pauli_state initial_state, double p){
        auto chpSim = ChpSimulator<TOTAL_QUBITS>(rand_gen);
        auto errorSim = FixedErrorSimPost<decltype(chpSim)>(chpSim, errorOp, error_n);
        auto idleSim = AutoIdlerSim<decltype(errorSim)>(errorSim);
        Planar::Planar<3, decltype(idleSim)> planar(errorSim, initial_state, p);

        // Planar::Planar<3, decltype(chpSim)> clone(planar.sim.get_error_free(), initial_state, p,  false);
        // {
        //     decltype(planar)::ErrorSyndrome s{};
        //     clone.perform_round_set(s);
        // }
        // planar.sim.get_error_free() = clone.sim;
        return planar;
}

int getSingleCount(){
    auto sim2 = getPlanarErrorSim(get_rand_gen(), single_qubit_pauli::I, -1, single_pauli_state::ZERO, 0.001);

    auto stabiliser4 = sim2.extractFTSyndrome();
    std::cout << sim2.sim.baseSim.getSingleCount() << " ";
    sim2.correctFTSyndrome(stabiliser4);
    return sim2.sim.baseSim.getSingleCount();
}

int getDoubleCount(){
    auto sim2 = getPlanarErrorSim(get_rand_gen(), single_qubit_pauli::I, -1, single_pauli_state::ZERO, 0.001);

    auto stabiliser4 = sim2.extractFTSyndrome();
    std::cout << sim2.sim.baseSim.getDoubleCount() << " ";
    // sim2.correctFTSyndrome(stabiliser4);
    return sim2.sim.baseSim.getDoubleCount();
}

std::vector<CodeTools::two_qubit_pauli> getTwoErrors() {
    std::vector<CodeTools::two_qubit_pauli> two_errors;
    for(single_qubit_pauli p1 : {single_qubit_pauli::I, single_qubit_pauli::X, single_qubit_pauli::Z, single_qubit_pauli::Y}){
        for(single_qubit_pauli p2 : {single_qubit_pauli::I, single_qubit_pauli::X, single_qubit_pauli::Z, single_qubit_pauli::Y}){
            if(!(p1 == single_qubit_pauli::I && p2 == single_qubit_pauli::I)) {
                two_errors.emplace_back(p1, p2);
            }
        }
    }
    return two_errors;
}

int main(){
    int single_count = getSingleCount();
    int double_count = getDoubleCount();
    std::cout << single_count << std::endl;
    std::cout << double_count << std::endl;

    auto twoErrors = getTwoErrors();

    auto rand_gen = get_rand_gen();

    int count = 0;
    int total_count = 0;
    for (int i = 0 ; i < single_count ; i++){
        for(auto error : {single_qubit_pauli::X,single_qubit_pauli::Y,single_qubit_pauli::Z}){
            for(int j = 0; j < 1000; j++){
                //Prepare the sim with injected errors.  
                auto sim2 = getPlanarErrorSim(rand_gen, error, i, single_pauli_state::PLUS, 0.001);

                auto stabiliser4 = sim2.extractFTSyndrome();
                sim2.correctFTSyndrome(stabiliser4);

                ChpSimulator<TOTAL_QUBITS> chp_reference(rand_gen);
                Planar::Planar<3, ChpSimulator<TOTAL_QUBITS>> sim3(chp_reference, single_pauli_state::PLUS, 0.0001, true);;

                auto stabiliser3 = sim3.perform_round();
                sim3.simple_correct(stabiliser3);

                if(sim3.sim.get_error_free() == sim2.sim.get_error_free()){
                    count++;
                } else {
                    sim2.sim.get_error_free().normalize();
                    std::cout << i << " " << error << " " << stabiliser4 << std::endl;
                }
                total_count ++;
            }
        }
    }

    for (int i = 0 ; i < double_count; i++){
        for(auto error : twoErrors){
            for(int j = 0; j < 1000; j++){
                //Prepare the sim with injected errors.  
                auto sim2 = getPlanarErrorSim(rand_gen, error, i, single_pauli_state::PLUS, 0.001);
                
                auto stabiliser4 = sim2.extractFTSyndrome();
                sim2.correctFTSyndrome(stabiliser4);

                ChpSimulator<TOTAL_QUBITS> chp_reference(rand_gen);
                Planar::Planar<3, ChpSimulator<TOTAL_QUBITS>> sim3(chp_reference, single_pauli_state::PLUS, 0.0001, true);;

                auto stabiliser3 = sim3.perform_round();
                sim3.simple_correct(stabiliser3);

                if(sim3.sim.get_error_free() == sim2.sim.get_error_free()){
                    count++;
                } else {
                    sim2.sim.get_error_free().normalize();
                    std::cout << i << " " << error.first  << " " << error.second << " " << stabiliser4 << std::endl;
                }
                total_count ++;
            }
        }
    }
    std::cout << count << " " << total_count << std::endl;
    
    return 0;
}