// Created by Alexis Shaw on 17/02/2021

#include <utility>
#include <array>
#include <vector>
#include <iostream>
// #include <algorithm>
// #include <map>
#include "sim/chp.h"
#include "BaseTypes.h"
#include "ErrorModels/TrivialErrorModel.h"
//#include "absl/container/flat_hash_map.h"
#include "PerfectMatching.h"
#include "boost/container/static_vector.hpp"
#include "boost/container/small_vector.hpp"
#include "prettyprint.h"

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>

namespace CodeTools::Planar {
    template<int distance>
    struct PlanarBase {
        public:
        static constexpr int LATTICE_DIMENSION{2*distance - 1};
        static constexpr int TOTAL_QUBITS{LATTICE_DIMENSION * LATTICE_DIMENSION};
        static constexpr int REQUIRED_ANCS{(TOTAL_QUBITS - 1)/2};
        static constexpr int REQUIRED_QUBITS{TOTAL_QUBITS - REQUIRED_ANCS};
        static constexpr int NUM_LOGICAL{1};

        using ErrorSyndrome = std::array<boost::container::static_vector<std::array<int,2>, REQUIRED_ANCS/2>,2>;
        using FTErrorSyndrome = std::array<boost::container::static_vector<std::array<int,3>, (REQUIRED_ANCS/2) * (distance + 1)>,2>;
        using InitialState = single_pauli_state; //Only X and Z supported for now.

        using ChangeArray = std::pair<std::array<std::array<int, REQUIRED_ANCS/2>, distance-1>,
                                      std::array<std::array<int, REQUIRED_ANCS/2>, distance-1>>;
        //Converts a coordinate to a qubit index
        inline static constexpr int get_qubit_index(int x, int y){
            return LATTICE_DIMENSION*x + y;
        }

        //gets the type of a qubit at a given qubit location
        inline static constexpr int get_type(int x, int y){
            return x%2;
        }
        inline static constexpr int get_type(int i){
            return (i/LATTICE_DIMENSION)%2;
        }

        inline static constexpr int get_type(std::pair<int, int>& location){
            return location.first%2;
        }

        inline static constexpr std::pair<int, int> get_qubit_location(int i){
            return std::make_pair(i / LATTICE_DIMENSION , i % LATTICE_DIMENSION);
        }

        static typename boost::container::small_vector<std::array<std::array<int,2>,2>, 10> match_planar_2D(ErrorSyndrome synd, int stab_type);
        static typename boost::container::small_vector<std::array<std::array<int,3>,2>, 10> match_planar_3D(FTErrorSyndrome synd, int stab_type);

        //returns a list of stabiliser indices, and stabiliser for a given type.
        static std::pair<boost::container::static_vector<int, 4>, short> get_stab(int i)
        {
            constexpr int row_length = 2 * distance - 1;
            int ancilla = 2*i+1;
            auto [row, column] = get_qubit_location(ancilla);
            int type = row % 2;
            using static_vector = boost::container::static_vector<int, 4>;

            if((row != 0) && (column != 0) && (row != (row_length - 1)) && (column != (row_length -1))){
                static_vector retval{ancilla - row_length, ancilla - 1,  ancilla + 1, ancilla + row_length};
                return std::make_pair(retval, type);
            }
            //top margin
            if (row == 0) {
                static_vector retval{-1, ancilla - 1, ancilla + 1, ancilla + row_length};
                return std::make_pair(retval, type);
            }
            //bottom margin
            if (row == row_length - 1) {
                static_vector retval{ancilla - row_length,  ancilla - 1, ancilla + 1, -1};
                return std::make_pair(retval, type);
            }
            //left margin
            if (column == 0) {
                static_vector retval{ancilla - row_length, -1, ancilla + 1, ancilla + row_length};
                return std::make_pair(retval, type);
            }
            //right margin
            if (column == row_length - 1) {
                static_vector retval{ancilla - row_length, ancilla - 1, -1, ancilla + row_length};
                return std::make_pair(retval, type);
            }
            //Control should never reach here
            return std::make_pair(static_vector(), -1);
        }
    };



    template<int distance, typename QubitSimulator=ChpSimulator<PlanarBase<distance>::TOTAL_QUBITS>>
    class Planar : public PlanarBase<distance>{
        public:
            using ErrorSyndrome = typename PlanarBase<distance>::ErrorSyndrome;
            using FTErrorSyndrome = typename PlanarBase<distance>::FTErrorSyndrome;
            using InitialState = typename PlanarBase<distance>::InitialState;
            using ChangeArray = typename PlanarBase<distance>::ChangeArray;

            using T_2DMatch = typename boost::container::small_vector<std::array<std::array<int,2>,2>, 10>;
            using T_3DMatch = typename boost::container::small_vector<std::array<std::array<int,3>,2>, 10>;

            static constexpr int LATTICE_DIMENSION{PlanarBase<distance>::LATTICE_DIMENSION};
            static constexpr int TOTAL_QUBITS{PlanarBase<distance>::TOTAL_QUBITS};
            static constexpr int REQUIRED_ANCS{PlanarBase<distance>::REQUIRED_ANCS};
            static constexpr int REQUIRED_QUBITS{PlanarBase<distance>::REQUIRED_QUBITS};
            static constexpr int NUM_LOGICAL{PlanarBase<distance>::NUM_LOGICAL};

            static std::pair<boost::container::static_vector<int, 4>, int> get_stab(int i){
                return PlanarBase<distance>::get_stab(i);
            }
            //gets the type of a qubit at a given qubit location
            inline static constexpr int get_type(int x, int y){
                return PlanarBase<distance>::get_type(x,y);
            }
            inline static constexpr int get_type(int i){
                return PlanarBase<distance>::get_type(i);
            }

            inline static constexpr int get_type(std::pair<int, int>& location){
                return PlanarBase<distance>::get_type(location);
            }

            //Converts a coordinate to a qubit index
            static constexpr int get_qubit_index(int x, int y){
                return PlanarBase<distance>::get_qubit_index(x,y);
            }

            static constexpr std::pair<int, int> get_qubit_location(int i){
                return PlanarBase<distance>::get_qubit_location(i);
            }
            static ChangeArray get_changes(FTErrorSyndrome &syndrome){
                return PlanarBase<distance>::get_changes(syndrome);
            }
            static T_2DMatch match_planar_2D(ErrorSyndrome synd, int stab_type){
                return PlanarBase<distance>::match_planar_2D(synd, stab_type);
            }
            static T_3DMatch match_planar_3D(FTErrorSyndrome synd, int stab_type){
                return PlanarBase<distance>::match_planar_3D(synd, stab_type);
            }

            static ErrorSyndrome getSyndromeFromInt(int measurement);

            void applyMatch2D(T_2DMatch&  matches, int stab_type);
            void applyMatch3D(T_3DMatch& matches, int stab_type);

            QubitSimulator sim{};
            const InitialState initialState;
        public:
            explicit Planar(const QubitSimulator& sim_in, const InitialState& state, bool do_init=true);
            explicit Planar(QubitSimulator&& sim_in, const InitialState& state, bool do_init=true);
            explicit Planar(QubitSimulator &sim_in, const InitialState& state,double alpha, double beta, double gamma, bool do_init=true);

        private:
            void doStabs();
        public:
            ErrorSyndrome perform_round();
            double perform_round_set(const ErrorSyndrome& i);

            void simple_correct(ErrorSyndrome error_syndrome){
                auto x_match = match_planar_2D(error_syndrome,0);
                auto z_match = match_planar_2D(error_syndrome, 1);

                applyMatch2D(x_match, 0);
                applyMatch2D(z_match, 1);
            }

            FTErrorSyndrome extractFTSyndrome();
            void correctFTSyndrome(FTErrorSyndrome error_syndrome){
                auto x_match = match_planar_3D(error_syndrome,0);
                auto z_match = match_planar_3D(error_syndrome, 1);

                applyMatch3D(x_match, 0);
                applyMatch3D(z_match, 1);
            }

            [[nodiscard]] bool confirm_no_error() const;

            explicit operator std::string();

            template<typename S>
            friend
            class FifteenSeven;
    };

    class weight_builder{
        public:
        struct graph_factory{
            const int width;
            const int height;
            const int thickness;

            int get_index(int i, int j, int k ){
                return i + (j*width) + k * height * width;
            }
            int num_nodes(){
                return width * height* thickness;
            }
            int get_x(int v){
            return v % width;
            }
            int get_y(int v){
            return (v/width)%height;
            }
            int get_z(int v){
            return v/(height * width);
            }
        };

        typedef adjacency_list < vecS, vecS, undirectedS, no_property, property<edge_weight_t, float>> graph_t;
        typedef graph_traits< graph_t >::vertex_descriptor vertex_descriptor;
        typedef std::pair<int, int> Edge;

        static graph_t get_graph(int width, int height, int thickness, bool s_type_is_x, double p){
            graph_factory fac {width,height,thickness};
            
            double p_a = 16.0* p *(1.0 / 15.0) * std::pow(1.0 - 4.0 * p * (1.0/ 15.0),3) * (1.0-p) +
                        p * std::pow(1.0 - 4.0 * p *(1.0/ 15.0),4);
            double p_b = 8.0 * p * (1.0 / 15.0) * (1.0 - 4.0 * p * (1.0/ 15.0)) * (1  - 2.0 * p * (1.0/3.0)) +
                        2.0 * p * (1.0/3.0) * std::pow(1.0 - 4.0 * p *(1.0/ 15.0),2);
            double p_c = 16.0* p *(1.0 / 15.0) * std::pow(1.0 - 4.0 * p * (1.0/ 15.0),3);
            double p_d = 16.0* p *(1.0 / 15.0) * std::pow(1.0 - 4.0 * p * (1.0/ 15.0),3) * std::pow(1 - 8.0*p*(1.0/15.0),2) * (1  - 2.0 * p * (1.0/3.0)) +
                        16.0* p *(1.0 / 15.0) * std::pow(1.0 - 4.0 * p * (1.0/ 15.0),4) * (1 - 8.0*p*(1.0/15.0)) * (1  - 2.0 * p * (1.0/3.0)) + 
                        2.0 * p * (1.0/3.0) * std::pow(1.0 - 4.0 * p * (1.0/ 15.0),4) * std::pow(1 - 8.0*p*(1.0/15.0),2);
            double p_e = 8.0 * p * (1.0 / 15.0) * (1.0 - 4.0 * p * (1.0/ 15.0));
            double p_f = p_e;
            float w_a = static_cast<float>(- std::log(p_a));
            float w_b = static_cast<float>(- std::log(p_b));
            float w_c = static_cast<float>(- std::log(p_c));
            float w_d = static_cast<float>(- std::log(p_d));
            float w_e = static_cast<float>(- std::log(p_e));
            float w_f = static_cast<float>(- std::log(p_f));

            const int num_nodes = fac.num_nodes();
            std::vector<Edge> edge_list{};
            std::vector<int> weights{};
            for(int i = 0; i < fac.width; i++){
            for(int j = 0; j < fac.height; j++){
                for(int k = 0; k < fac.thickness; k++){
                //Downward Edge - A
                if(k > 0){
                    edge_list.push_back(Edge(fac.get_index(i,j,k), fac.get_index(i,j,k - 1)));
                    // std::cout << i << " "<< j << " "<< k << "--" << i << " "<< j<< " "<< k-1<< "\n";
                    // std::cout << fac.get_index(i,j,k) << "--"<< fac.get_index(i-1,j,k) << "\n";
        //               weights.push_back(1.0f);
                    weights.push_back(w_a);
                }
                //Northward Edge - B
                if(j > 0){
                    edge_list.push_back(Edge(fac.get_index(i,j,k), fac.get_index(i,j - 1,k)));
                    // std::cout << i << " "<< j << " "<< k << "--" << i << " "<< j-1<< " "<< k <<"\n";
                    // std::cout << fac.get_index(i,j,k) << "--"<< fac.get_index(i-1,j,k) << "\n";
        //               weights.push_back(1.0f);
                    weights.push_back(w_b);
                }
                //North-Downward Edge - C
                if(k > 0 && j > 0){
                    edge_list.push_back(Edge(fac.get_index(i,j,k), fac.get_index(i,j - 1,k-1)));
                    // std::cout << i << " "<< j << " "<< k << "--" << i << " "<< j-1<< " "<< k-1 <<"\n";
                    // std::cout << fac.get_index(i,j,k) << "--"<< fac.get_index(i-1,j,k-1) << "\n";
        //               weights.push_back(1.0f);
                    weights.push_back(w_c);
                }
                //Westward Edge -- D
                if(i > 0){
                    edge_list.push_back(Edge(fac.get_index(i,j,k), fac.get_index(i - 1,j,k)));
                    // std::cout << i << " "<< j << " "<< k << "--" << i-1 << " "<< j<< " "<< k << "\n";
                    // std::cout << fac.get_index(i,j,k) << "--"<< fac.get_index(i-1,j,k) << "\n";
        //               weights.push_back(1.0f);
                    weights.push_back(w_d);
                }
                //West-Downward Edge -- E
                if(i > 0 && k > 0){
                    edge_list.push_back(Edge(fac.get_index(i,j,k), fac.get_index(i - 1,j,k-1)));
                    // std::cout << i << " "<< j << " "<< k << "--" << i-1 << " "<< j<< " "<< k-1 << "\n";
                    // std::cout << fac.get_index(i,j,k) << "--"<< fac.get_index(i-1,j,k-1) << "\n";
        //               weights.push_back(1.0f);
                    weights.push_back(w_e);
                }
                //North-Down-Eastward Edge
                if(i < (fac.width - 1) && j > 0 && k > 0){
                    edge_list.push_back(Edge(fac.get_index(i,j,k), fac.get_index(i + 1,j-1,k-1)));
                    //std::cout << i << " "<< j << " "<< k << "--" << i+1 << " "<< j-1<< " "<< k-1 << "\n";
                    // std::cout << fac.get_index(i,j,k) << "--"<< fac.get_index(i+1,j-1,k-1) << "\n";
        //               weights.push_back(1.0f);
                    weights.push_back(w_f);
                }
                }
            }
            }
            if(s_type_is_x){
                for(int i = 0; i < fac.width; i++){
                    for(int k = 0; k < fac.thickness; k++){
                        //Northward Boundary
                        edge_list.push_back(Edge(fac.get_index(i,0,k), fac.num_nodes()));
        //                 weights.push_back(1.0f);
                        weights.push_back(w_b);
                        //Southern Boundary
                        edge_list.push_back(Edge(fac.get_index(i,fac.height - 1,k), fac.num_nodes() + 1));
        //                 weights.push_back(1.0f);
                        weights.push_back(w_b);
                    }
                }
            } else {
                for(int j = 0; j < fac.height; j++){
                    for(int k = 0; k < fac.thickness; k++){
                        //Westward Boundary
                        edge_list.push_back(Edge(fac.get_index(0,j,k), fac.num_nodes()));
        //                 weights.push_back(1.0f);
                        weights.push_back(w_d);
                        //Eastward Boundary
                        edge_list.push_back(Edge(fac.get_index(fac.width - 1,j,k), fac.num_nodes() + 1));
        //                 weights.push_back(1.0f);
                        weights.push_back(w_d);
                    }
                }
            }
            
            
            graph_t g(edge_list.begin(), edge_list.end(), weights.begin(), num_nodes + 2);
            return g;
        }
        
        static std::pair<std::vector<std::vector<float>>, std::vector<std::vector<vertex_descriptor>>> get_distances(graph_t& g){
            std::vector<std::vector<float>> d_l{};
            std::vector<std::vector<vertex_descriptor>> p_l{};
            for (vertex_descriptor i = 0; i < num_vertices(g); ++i){
            std::vector<vertex_descriptor> p(num_vertices(g));
            std::vector<float> d(num_vertices(g));

            vertex_descriptor s = vertex(i, g);
            dijkstra_shortest_paths(g, s,
                                    predecessor_map(boost::make_iterator_property_map(p.begin(), get(boost::vertex_index, g))).
                                    distance_map(boost::make_iterator_property_map(d.begin(), get(boost::vertex_index, g))));

            d_l.push_back(std::move(d));
            p_l.push_back(std::move(p));
            }
            return std::make_pair(d_l, p_l);
        }

    };

    template<int distance>
    typename boost::container::small_vector<std::array<std::array<int,2>,2>, 10>
    PlanarBase<distance>::match_planar_2D(PlanarBase::ErrorSyndrome synd, int stab_type) {
        if(synd[stab_type].size() == 0)
            return {};
        size_t N = synd[stab_type].size();
        PerfectMatching pm(2*N, 4*N*N);
        pm.options.verbose = false;

        std::sort(synd[stab_type].begin(), synd[stab_type].end(), [](std::array<int,2> a, std::array<int,2> b){
            if(a[0] != b[0]){
                return a[0] > b[0];
            } else {
                return a[1] > b[1];
            }
        });

        boost::container::small_vector<std::array<int,2>,10> boundary_node_positions{};
        //std::cout << N << std::endl;
        //Add complete graph of syndromes;
        for(size_t i = 0; i < N; ++i) {
            auto[p0, p1] = synd[stab_type][i];
            for (size_t j = i + 1; j < N; ++j) {
                auto[q0, q1] = synd[stab_type][j];
                //std::cout << i << " " << j << " " << abs(p0 - q0) + abs(p1 - q1) << std::endl;
                pm.AddEdge(i,j,abs(p0 - q0) + abs(p1 - q1));
                //std::cout << i+N << " " << j+N << " " << 0 << std::endl;
                pm.AddEdge(i+N,j+N,0);
            }
        }
        for(size_t i = 0; i < N; ++i) {
            auto[p0, p1] = synd[stab_type][i];
            int r0 = p0, r1 = p1;
            if(stab_type == 1){
                r0 = (p0 < distance)? -1: LATTICE_DIMENSION;
                //std::cout << i << " " << i + N << " " << abs(p0 - r0)<< std::endl;
                pm.AddEdge(i, i + N, abs(p0 - r0));
            } else {
                r1 = (p1 < distance)? -1: LATTICE_DIMENSION;
                //std::cout << i << " " << i + N << " " << abs(p1 - r1)<< std::endl;
                pm.AddEdge(i, i + N, abs(p1 - r1));
            }
            boundary_node_positions.push_back({r0,r1});
        }

        pm.Solve();
        //std::cout << "Matching results" << std::endl;
        boost::container::small_vector<std::array<std::array<int,2>,2>, 10> matching{};
        for(size_t i = 0; i< N ; ++i){
            size_t j = pm.GetMatch(i);
            if(j < i){
                //Only do one copy of each matching;
                matching.push_back({synd[stab_type][i],synd[stab_type][j]}); // Then add that matching to the list
            } else if(j >= N){
                //Boundary node
                matching.push_back({synd[stab_type][i],boundary_node_positions[j - N]});
            }
            //std::cout << i << " " << j << std::endl;
        }
        return matching;
    }

    template<int distance>
    typename boost::container::small_vector<std::array<std::array<int,3>,2>, 10>
    PlanarBase<distance>::match_planar_3D(PlanarBase::FTErrorSyndrome synd, int stab_type) {
        if(synd[stab_type].size() == 0)
            return {};
        size_t N = synd[stab_type].size();
        PerfectMatching pm(2*N, 4*N*N);
        pm.options.verbose = false;

        std::sort(synd[stab_type].begin(), synd[stab_type].end(), [](std::array<int,3> a, std::array<int,3> b){
            if(a[0] != b[0]){
                return a[0] > b[0];
            } else if (a[1] != b[1]) {
                return a[1] > b[1];
            } else {
                return a[2] > b[2];
            }
        });

        boost::container::small_vector<std::array<int,3>, 10> boundary_node_positions{};
        //Add complete graph of syndromes;
        for(size_t i = 0; i < N; ++i) {
            auto[p0, p1, pt] = synd[stab_type][i];
            for (size_t j = i + 1; j < N; ++j) {
                auto[q0, q1, qt] = synd[stab_type][j];
                pm.AddEdge(i,j,abs(p0 - q0) + abs(p1 - q1) + abs(pt - qt));
                pm.AddEdge(i+N,j+N,0);
            }
        }
        for(size_t i = 0; i < N; ++i) {
            auto [p0, p1, pt] = synd[stab_type][i];
            int r0 = p0, r1 = p1, rt = pt;
            if(stab_type == 1){
                r0 = (p0 < distance)? -1: LATTICE_DIMENSION;
                //std::cout << i << " " << i + N << " " << abs(p0 - r0)<< std::endl;
                pm.AddEdge(i, i + N, abs(p0 - r0));
            } else {
                r1 = (p1 < distance)? -1: LATTICE_DIMENSION;
                //std::cout << i << " " << i + N << " " << abs(p1 - r1)<< std::endl;
                pm.AddEdge(i, i + N, abs(p1 - r1));
            }
            boundary_node_positions.push_back({r0,r1, rt});
        }

        pm.Solve();
        //std::cout << "Matching results" << std::endl;
        boost::container::small_vector<std::array<std::array<int,3>,2>, 10> matching{};
        for(size_t i = 0; i< N ; ++i){
            size_t j = pm.GetMatch(i);
            if(j < i){
                //Only do one copy of each matching;
                matching.push_back({synd[stab_type][i],synd[stab_type][j]}); // Then add that matching to the list
            } else if(j >= N){
                //Boundary node
                matching.push_back({synd[stab_type][i],boundary_node_positions[j - N]});
            }
            //std::cout << i << " " << j << std::endl;
        }
        return matching;
    }

    template <int distance, typename QubitSimulator>
    inline Planar<distance, QubitSimulator>::Planar(const QubitSimulator& sim_in, const InitialState& state, bool do_init):
    sim{sim_in}, initialState{state}{
        sim.enable_errors = 0;
        if(do_init){
            switch(state){
                case single_pauli_state::ZERO:
                    break;
                case single_pauli_state::ONE:
                    for(int i = 0; i <= distance; i++){
                        sim.x(get_qubit_index(i,i));
                    }
                    break;
                case single_pauli_state::PLUS:
                    for(int i = 0; i < TOTAL_QUBITS; i++){
                        if((i%2) == 0)
                            sim.h(i);
                    }
                    break;
                case single_pauli_state::MINUS:
                    for(int i = 0; i < TOTAL_QUBITS; i++){
                        if((i%2) == 0) {
                            sim.h(i);
                        }
                    }
                    for(int i = 0; i <= distance; i++){
                        sim.z(get_qubit_index(i,i));
                    }
                    break;
                case single_pauli_state::I_MINUS:
                    assert(false);// not supported for Surface code initialisation.
                    break;
                case single_pauli_state::I_PLUS:
                    assert(false);// not supported for Surface code initialisation.
                    break;
            }
        }
        sim.enable_errors = 1;
        //std::cout << "constructor 1" << std::endl;

    }

    template <int distance, typename QubitSimulator>
    inline Planar<distance, QubitSimulator>::Planar(QubitSimulator &&sim_in, const InitialState& state, bool do_init):
    sim{std::move(sim_in)}, initialState{state}{
        //std::cout << "constructor2" << std::endl;
        sim.enable_errors = 0;
        if(do_init){
            switch(state){
                case single_pauli_state::ZERO:
                    break;
                case single_pauli_state::ONE:
                    for(int i = 0; i <= distance; i++){
                        sim.x(get_qubit_index(i,i));
                    }
                    break;
                case single_pauli_state::PLUS:
                    for(int i = 0; i < TOTAL_QUBITS; i++){
                        if((i%2) == 0)
                            sim.h(i);
                    }
                    break;
                case single_pauli_state::MINUS:
                    for(int i = 0; i < TOTAL_QUBITS; i++){
                        if((i%2) == 0) {
                            sim.h(i);
                        }
                    }
                    for(int i = 0; i <= distance; i++){
                        sim.z(get_qubit_index(i,i));
                    }
                    break;
                case single_pauli_state::I_MINUS:
                    assert(false);// not supported for Surface code initialisation.
                    break;
                case single_pauli_state::I_PLUS:
                    assert(false);// not supported for Surface code initialisation.
                    break;
            }
        }
        sim.enable_errors = 1;
    }

    template <int distance, typename QubitSimulator>
    inline Planar<distance, QubitSimulator>::Planar(QubitSimulator &sim_in, const InitialState& state,double alpha, double beta, double gamma, bool do_init):
            sim{std::move(sim_in)}, initialState{state}{
        //std::cout << "constructor2" << std::endl;
        sim.enable_errors = 0;
        if(do_init){
            switch(state){
                case single_pauli_state::ZERO:
                    break;
                case single_pauli_state::ONE:
                    for(int i = 0; i <= distance; i++){
                        sim.x(get_qubit_index(i,i));
                    }
                    break;
                case single_pauli_state::PLUS:
                    for(int i = 0; i < TOTAL_QUBITS; i++){
                        if((i%2) == 0)
                            sim.h(i);
                    }
                    break;
                case single_pauli_state::MINUS:
                    for(int i = 0; i < TOTAL_QUBITS; i++){
                        if((i%2) == 0) {
                            sim.h(i);
                        }
                    }
                    for(int i = 0; i <= distance; i++){
                        sim.z(get_qubit_index(i,i));
                    }
                    break;
                case single_pauli_state::I_MINUS:
                    assert(false);// not supported for Surface code initialisation.
                    break;
                case single_pauli_state::I_PLUS:
                    assert(false);// not supported for Surface code initialisation.
                    break;
            }
        }
        for(int i = 0; i < TOTAL_QUBITS; ++i) {
            if((i%2) == 0){
                sim.y(i, alpha);
                sim.z(i, beta);
                sim.y(i, gamma);
            }
        }
        sim.enable_errors = 1;
    }

    template<int distance, typename QubitSimulator>
    inline void Planar<distance, QubitSimulator>::doStabs(){
        //X-Stabilisers
        std::array<decltype(get_stab(0)), REQUIRED_ANCS> stabs{};
        for (int j=0; j<REQUIRED_ANCS; j++){
            stabs[j] = get_stab(j);
        }
        for (int j=0; j<REQUIRED_ANCS; j++) {
            int ancilla = 2 * j + 1;
            if(stabs[j].second == 0)
                sim.h(ancilla);
        }
        for (int i = 0; i < 4 ; i++){
            for (int j=0; j<REQUIRED_ANCS; j++) {
                auto &[stab, type] = stabs[j];
                int ancilla = 2 * j + 1;
                if (type == 1 && (stab[i] != -1)) {
                    sim.cx(stab[i], ancilla);
                }
            }
        }
        for (int i = 0; i < 4 ; i++){
            for (int j=0; j<REQUIRED_ANCS; j++) {
                auto &[stab, type] = stabs[j];
                int ancilla = 2 * j + 1;
                if (type == 0 && (stab[i] !=-1)) {
                    sim.cx(ancilla, stab[i]);
                }
            }
        }

        for (int j=0; j<REQUIRED_ANCS; j++) {
            int ancilla = 2 * j + 1;
            if(stabs[j].second == 0)
                sim.h(ancilla);
        }
    }

    template <int distance, typename QubitSimulator>
    inline typename Planar<distance, QubitSimulator>::ErrorSyndrome Planar<distance, QubitSimulator>::getSyndromeFromInt(int measurement){
        boost::container::small_vector<measurement_state, 20> retval{};
        for (size_t i = 0; i < REQUIRED_ANCS; ++i){
            retval.push_back((measurement & (1<<i)? measurement_state::RANDOM_TRUE : measurement_state::RANDOM_FALSE));
        }
        auto measurement_results = retval;

        boost::container::small_vector<int,10> ancilla_nos{};
        for (int i = 0; i < REQUIRED_ANCS; ++i) {
            ancilla_nos.push_back(2 * i + 1);
        }

        ErrorSyndrome syndrome{};
        for (int i = 0; i < REQUIRED_ANCS; i++) {
            auto[stabs, type] = get_stab(i);
            if (measurement_results[i] == measurement_state::RANDOM_TRUE) {
                auto[x, y] = get_qubit_location(ancilla_nos[i]);
                syndrome[type].push_back({x, y});
            }
        }
        return syndrome;
    }

    template <int distance, typename QubitSimulator>
    inline typename Planar<distance, QubitSimulator>::ErrorSyndrome Planar<distance, QubitSimulator>::perform_round() {

        //X-Stabilisers
        doStabs();

        boost::container::small_vector<int,10> ancilla_nos{};
        for (int i = 0; i < REQUIRED_ANCS; ++i) {
            ancilla_nos.push_back(2 * i + 1);
        }

        auto measurement_results = sim.measure_all(ancilla_nos);

        ErrorSyndrome syndrome{};
        for (int i = 0; i < REQUIRED_ANCS; i++) {
            auto[stabs, type] = get_stab(i);
            if (measurement_results[i] == measurement_state::RANDOM_TRUE) {
                auto[x, y] = get_qubit_location(ancilla_nos[i]);
                syndrome[type].push_back({x, y});
            }
        }

        for (int i = 0; i < 2; ++i) {
            for (size_t j = 0; j < syndrome[i].size(); ++j) {
                sim.x(get_qubit_index(syndrome[i][j][0], syndrome[i][j][1])); //reset all 1-qubits.
            }
        }
        return syndrome;
    }

    template <int distance, typename QubitSimulator>
    inline double Planar<distance, QubitSimulator>::perform_round_set(const ErrorSyndrome& syndrome){
        //Do stab circuits;
        doStabs();

        boost::container::small_vector<int, 10> ancilla_nos{};
        for(int i = 0; i<REQUIRED_ANCS; ++i){
            ancilla_nos.push_back(2*i + 1);
        }

        boost::container::small_vector<int, 10> true_qubits{};
        for(const auto& type_syndromes : syndrome) for(const auto& loc : type_syndromes)
            true_qubits.push_back(get_qubit_index(loc[0],loc[1]));

        std::sort(true_qubits.begin(), true_qubits.end());

        boost::container::small_vector<measurement_state, 10> measurements;

        for(int ancilla_no : ancilla_nos){
            if (std::find(true_qubits.begin(), true_qubits.end(), ancilla_no) == true_qubits.end()){
                measurements.push_back(measurement_state::RANDOM_FALSE);
            } else {
                measurements.push_back(measurement_state::RANDOM_TRUE);
            }
        }

        double retval = sim.measure_all_set(ancilla_nos, measurements);

        //Reset
        for(int i=0; i < 2; ++i){
            for(size_t j = 0 ; j < syndrome[i].size(); ++j){
                sim.x(get_qubit_index(syndrome[i][j][0], syndrome[i][j][1])); //reset all 1-qubits.
            }
        }
        //std::cout << sim.reg.steps[1] << std::flush;
        return retval;
    }

    template <int distance, typename QubitSimulator>
    inline typename Planar<distance, QubitSimulator>::FTErrorSyndrome Planar<distance, QubitSimulator>::extractFTSyndrome(){
        FTErrorSyndrome FTSyndrome{};
        ErrorSyndrome prev{};
        for(int t = 0; t < distance; ++t){
            //std::cout << "t" << t << std::flush;
            std::array<std::array<std::array<bool, LATTICE_DIMENSION>, LATTICE_DIMENSION>,2> parity{};
            ErrorSyndrome round = perform_round();
            //std::cout << std::endl;
//            std::cout << round << std::endl;
            for(int stab =0; stab < 2; stab++) {
                for (const auto&[x, y] : round[stab]) {
                    parity[stab][x][y] ^= 1;
                }
                for (const auto&[x, y] : prev[stab]) {
                    parity[stab][x][y] ^= 1;
                }
                for(int i = 0; i < LATTICE_DIMENSION ; ++i){
                    for(int j = 0; j < LATTICE_DIMENSION; ++j){
                        if(parity[stab][i][j]){
                            FTSyndrome[stab].push_back({i,j,t});
                        }
                    }
                }
            }
            prev = round;
        }
        sim.enable_errors = 0;
        {
            int t = distance;
            std::array<std::array<std::array<bool, LATTICE_DIMENSION>, LATTICE_DIMENSION>,2> parity{};
            ErrorSyndrome round = perform_round();
            //std::cout << round << std::endl;
            for(int stab =0; stab < 2; stab++) {
                for (const auto&[x, y] : round[stab]) {
                    parity[stab][x][y] ^= 1;
                }
                for (const auto&[x, y] : prev[stab]) {
                    parity[stab][x][y] ^= 1;
                }

                for(int i = 0; i < LATTICE_DIMENSION ; ++i){
                    for(int j = 0; j < LATTICE_DIMENSION; ++j){
                        if(parity[stab][i][j]){
                            FTSyndrome[stab].push_back({i,j,t});
                        }
                    }
                }
            }
            prev = round;
        }
        sim.enable_errors = 1;
        return FTSyndrome;
    }

    template <int distance, typename QubitSimulator>
    inline bool Planar<distance, QubitSimulator>::confirm_no_error() const{
        return false; // Not Implemented.

    }

    template <int distance, typename QubitSimulator>
    void Planar<distance, QubitSimulator>::applyMatch2D(boost::container::small_vector<std::array<std::array<int,2>,2>, 10>&  matches, int stab_type){
        boost::container::small_vector<int, 10> flips{};
        for (const auto m : matches) {
            const auto &[p0, p1] = m[0];
            const auto &[q0, q1] = m[1];

            for (int i = std::min(p0, q0); i < std::max(p0, q0) + 1; ++i){
                if(((i + p1) % 2) == 0) {
                    flips.push_back(get_qubit_index(i, p1));
                }
            }
            for (int i = std::min(p1, q1); i < std::max(p1, q1) + 1; ++i){

                if(((i + q0) % 2) == 0) {
                    flips.push_back(get_qubit_index(q0, i));
                }
            }
        }
        for(int q : flips){
            if(stab_type == 0){
                sim.z(q);
            } else {
                sim.x(q);
            }
        }
    }

    template <int distance, typename QubitSimulator>
    inline Planar<distance, QubitSimulator>::operator std::string(){
        return std::string("String Conversion Not Implemented yet");
    }

    template<int distance, typename QubitSimulator>
    void Planar<distance, QubitSimulator>::applyMatch3D(boost::container::small_vector<std::array<std::array<int,3>,2>, 10> &matches,
                                                        int stab_type) {
        boost::container::small_vector<int, 10> flips{};
        for (const auto m : matches) {
            const auto &[p0, p1, pt] = m[0];
            const auto &[q0, q1, qt] = m[1];

            for (int i = std::min(p0, q0); i < std::max(p0, q0) + 1; ++i){
                if(((i + p1) % 2) == 0) {
                    flips.push_back(get_qubit_index(i, p1));
                }
            }
            for (int i = std::min(p1, q1); i < std::max(p1, q1) + 1; ++i){

                if(((i + q0) % 2) == 0) {
                    flips.push_back(get_qubit_index(q0, i));
                }
            }
        }
        std::array<bool,TOTAL_QUBITS> flip_array{};
        for(int q : flips){
            flip_array[q] ^= 1;
        }
        for(size_t i = 0 ; i < flip_array.size(); ++i){
            if(flip_array[i]) {
                if (stab_type == 0) {
                    //std::cout << "c_z" << i <<std::flush;
                    sim.z(i);
                } else {
                    //std::cout << "c_x" << i <<std::flush;
                    sim.x(i);
                }
            }
        }
        //std::cout << std::endl;

    }
}
