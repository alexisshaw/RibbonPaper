#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <vector>
#include <array>
#include "PerfectMatching.h"
#include "boost/container/small_vector.hpp"
#include <boost/config.hpp>
#include <iostream>
#include <fstream>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>
#include <tuple>

// ----------------
// Function that does a correction in python
// ----------------
namespace py = pybind11;
using namespace boost;

using FTErrorSyndrome = std::vector<std::array<short, 3>>;
using FTMatching = boost::container::small_vector<std::array<std::array<short,3>,2>, 1024>;

float lookup_distance(py::array_t<float, py::array::c_style |py::array::forcecast>& distances, 
                       bool stab_type_x, int p0, int p1, int pt, int q0, int q1, int qt){
    auto r = distances.unchecked<6>();
    if(stab_type_x){
        return r(pt, p0/2, p1/2, qt, (q0 + 2)/2, q1/2);
    } else {
        return r(pt, p0/2, p1/2, qt, q0/2,(q1 + 2)/2);
    }
    
}

FTMatching match_planar_3D(const FTErrorSyndrome& synd, const bool stab_type_x, const short dx, const short dz, py::array_t<float, py::array::c_style |py::array::forcecast>& distances) {
    if(synd.size() == 0)
        return {};
    size_t N = synd.size();
    PerfectMatching pm(2*N, N*N + 2*N);
    pm.options.verbose = false;
    // py::print("Matching");

    boost::container::small_vector<std::array<short,3>, 10> boundary_node_positions{};
    //Add complete graph of syndromes;
    for(size_t i = 0; i < N; ++i) {
        auto[pt, p0, p1] = synd[i];
        for (size_t j = i + 1; j < N; ++j) {
            auto[qt, q0, q1] = synd[j];
            pm.AddEdge(i,j,lookup_distance(distances, stab_type_x, p0, p1, pt, q0, q1, qt));
            pm.AddEdge(i+N,j+N,0);
        }
    }

    for(size_t i = 0; i < N; ++i) {
        auto [pt, p0, p1] = synd[i];
        short r0 = p0, r1 = p1, rt = pt;
        if(stab_type_x){
            int d0 = lookup_distance(distances, stab_type_x, p0, p1, pt, -1, 0, 0);
            int dn = lookup_distance(distances, stab_type_x, p0, p1, pt, (2*dz - 1), 0, 0);
            r0 = (d0 < dn)? -1: (2*dz - 1);
            pm.AddEdge(i,i+N,std::min(d0, dn));
        } else {
            int d0 = lookup_distance(distances, stab_type_x, p0, p1, pt, 0, -1, 0);
            int dn = lookup_distance(distances, stab_type_x, p0, p1, pt, 0, (2*dx - 1), 0);
            r1 = (d0 < dn)? -1: (2*dx - 1);
            pm.AddEdge(i,i+N, std::min(d0, dn));
        }
        boundary_node_positions.push_back({rt,r0, r1});
    }

    pm.Solve();
    //std::cout << "Matching results" << std::endl;
    FTMatching matching{};
    for(size_t i = 0; i< N ; ++i){
        size_t j = pm.GetMatch(i);
        if(j < i){
            //Only do one copy of each matching;
            matching.push_back({synd[i],synd[j]}); // Then add that matching to the list
        } else if(j >= N){
            //Boundary node
            matching.push_back({synd[i],boundary_node_positions[j - N]});
        }
        //std::cout << i << " " << j << std::endl;
    }
    return matching;
}

int get_qubit_index(short p0, short p1, short dx, short dz){
    return ((2*dx - 1) * p0 + p1)/2;
    // return ((2*dx - 1) * p0 + p1);
}

void applyMatch3D(const FTMatching &matches, uint8_t* flip_array, short dx, short dz) {
    for (const auto m : matches) {
        const auto &[pt, p0, p1] = m[0];
        const auto &[qt, q0, q1] = m[1];

        for (int i = std::min(p0, q0); i < std::max(p0, q0) + 1; ++i){
            if(((i + p1) % 2) == 0) {
                flip_array[get_qubit_index(i, p1, dx, dz)] ^= 1;
            }
        }
        for (int i = std::min(p1, q1); i < std::max(p1, q1) + 1; ++i){
            if(((i + q0) % 2) == 0) {
                flip_array[get_qubit_index(q0, i, dx, dz)] ^= 1;
            }
        }
    }
}

void do_correction(const FTErrorSyndrome& axions, short dx, short dz, bool stab_type_x, uint8_t* pauli_frame, py::array_t<float, py::array::c_style |py::array::forcecast>& distances){
    auto matches = match_planar_3D(axions, stab_type_x, dx, dz, distances);

    // for(auto pair : matches){
    //     py::print("(", pair[0][0],", ", pair[0][1],", ", pair[0][2], ") <-> (", pair[1][0],", ", pair[1][1],", ", pair[1][2],")");
    // }

    applyMatch3D(matches, pauli_frame, dx, dz);
    return;
}

py::array get_corrections(py::array_t<short, py::array::c_style | py::array::forcecast> array, int64_t num_samples, short dx, short dz, short n_repeats, std::string stab_type, py::array_t<float, py::array::c_style |py::array::forcecast> distances){
    // ssize_t indim = array.ndim();
    // py::print("ndim",indim);
    // py::print("shape");
    // for (int i = 0; i < indim; ++i){
    //     py::print("d_",i, " ",array.shape()[i]);
    // }
    // py::print("span");
    // for (int i = 0; i < indim; ++i){
    //     py::print("d_",i, " ",array.strides()[i]);
    // }
    // py::print("size", array.size());

    bool s_type_is_x = (stab_type == "X" || stab_type =="x")? true:false;
    // py::print("stab_type = X", s_type_is_x);


    if(array.ndim() != 2)
        throw std::runtime_error("Input should be a 2-D NumPy array");
    if(array.shape()[1] != 4)
        throw std::runtime_error("Input shouls have size [N,4]");
    
    if(distances.ndim() != 6)
        throw std::runtime_error("distances map should be a 6-D NumPy array");
//     if(distances.shape()[0] != distances.shape()[3] ||
//        distances.shape()[1] != distances.shape()[4] ||
//        distances.shape()[2] != distances.shape()[5])
//         throw std::runtime_error("distances map should map between two 3-D points");
//     if(distances.shape()[0] != n_repeats)
//         throw std::runtime_error("distances map time-length does not match given time_span");
//     if((s_type_is_x    && (distances.shape()[2] != dx || distances.shape()[1] != (dz-1)))||
//        ((!s_type_is_x) && (distances.shape()[2] != (dx-1) || distances.shape()[1] != dz)))
//         throw std::runtime_error("distances map dimensions wrong for syndrome type and code distances");

    ssize_t stride[2] = {array.strides()[0], array.strides()[1]};
    auto* data = array.data();
    // for(int i = 0; i < array.shape()[0]; i++){
    //     const int64_t* inner_data = data + i * stride[0]/sizeof(int64_t);
    //     py::print("Row, ", i , ", ", inner_data[0], " ", inner_data[1], " ", inner_data[2], " ", inner_data[3]);
    // }

    int num_logical = dx*dz + (dx - 1)*(dz-1);
    std::vector<std::pair<const short*,const short*>> indicies(num_samples);
    // int num_logical = (2*dx - 1) * (2*dz - 1);
    for(auto [i,inner_data] = std::tuple{ssize_t{0}, array.data()};
        i < num_samples;
        i++){
        // Insert the axions for the given sample into a local vector;
        indicies[i].first = inner_data;
        for( ; (inner_data[0] <= i) && (inner_data < (data + array.size())); inner_data += stride[0]/sizeof(short)){
        }
        indicies[i].second = inner_data;
    }
    
    std::vector<uint8_t> result_vec(num_samples * num_logical);
    {   
        ssize_t s0 = stride[0];
        FTErrorSyndrome inner_axions{};
        //#pragma omp parallel for default(none) shared(indicies, result_vec) firstprivate(s0, inner_axions, dx, dz, s_type_is_x, num_samples, num_logical) schedule(dynamic, 50)
        for(int i=0 ; i < num_samples; ++i){
            // Insert the axions for the given sample into a local vector;
            for(const short* inner_data = indicies[i].first; inner_data < indicies[i].second; inner_data += s0/sizeof(short)){
                if(s_type_is_x){
                    inner_axions.push_back({{inner_data[1], short(inner_data[2] * 2 + 1), short(inner_data[3] * 2)}});
                } else {
                    inner_axions.push_back({{inner_data[1], short(inner_data[2] * 2), short(inner_data[3] * 2 + 1)}});
                }
            }
            // Print the axiom vector for this sample
            // py::print(i);
            // for(const auto& row: inner_axions){
            //     py::print(row[0], " ", row[1], " ", row[2]);
            // }

            do_correction(inner_axions, dx, dz, s_type_is_x, &(result_vec[i*num_logical]), distances);
            inner_axions.clear();
        }
    }

    // call pure C++ function

    // ssize_t ndim = 3;
    // std::vector<ssize_t> shape  = {result_vec.size()/num_logical, (2*dz - 1), (2*dx - 1)};
    // std::vector<ssize_t> strides = {sizeof(uint8_t)* num_logical,sizeof(uint8_t) * (2*dx - 1),  sizeof(uint8_t)};

    ssize_t ndim = 2;
    std::vector<ssize_t> shape  = {static_cast<ssize_t>(result_vec.size())/num_logical, num_logical};
    std::vector<ssize_t> strides = {static_cast<ssize_t>(sizeof(uint8_t))* num_logical, static_cast<ssize_t>(sizeof(uint8_t))};

    return py::array(py::buffer_info(
        result_vec.data(), //data as contiguous array
        sizeof(uint8_t), // size of a scalar
        py::format_descriptor<uint8_t>::format(), //data type
        ndim, shape, strides
    ));
}

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

py::array get_distances(int width, int height, int thickness, std::string stab_type, double p){
    bool s_type_is_x = (stab_type == "X" || stab_type =="x")? true:false;
    
    weight_builder::graph_t g = weight_builder::get_graph(width,height,thickness,s_type_is_x, p);

    auto[d_l, p_l] = weight_builder::get_distances(g);

    std::vector<ssize_t> shape  = {thickness, height, width, thickness, height, width};
    
    if(s_type_is_x){
        shape[4] += 2;
    } else {
        shape[5] += 2;
    }
    
    auto array = py::array_t<float>(shape);
    
    auto r = array.mutable_unchecked<6>();
    
    weight_builder::graph_factory fac {width,height,thickness};
    
    for(int i0 = 0; i0 < thickness; i0++){
        for(int j0 = 0; j0 < height; j0++){
            for(int k0 = 0; k0 < width; k0++){
                for(int i1 = 0; i1< thickness; i1++){
                    for(int j1 = 0; j1 < height; j1++){
                        for(int k1 = 0; k1 < width; k1++){
                            if(s_type_is_x){
                                r(i0, j0, k0, i1, 0, k1) = 0.0f;
                                r(i0, j0, k0, i1, height+1, k1) = 0.0f;
                                r(i0, j0, k0, i1, j1 + 1, k1) = d_l[fac.get_index(k0,j0,i0)][fac.get_index(k1,j1,i1)];
                            } else {
                                r(i0, j0, k0, i1, j1, 0) = 0.0f;
                                r(i0, j0, k0, i1, j1, width+1) = 0.0f;
                                r(i0, j0, k0, i1, j1, k1 + 1) = d_l[fac.get_index(k0,j0,i0)][fac.get_index(k1,j1,i1)];
                            }
                        }
                    }
                }
            }
        }
    }
    
    for(int i0 = 0; i0 < thickness; i0++){
        for(int j0 = 0; j0 < height; j0++){
            for(int k0 = 0; k0 < width; k0++){
                if(s_type_is_x){
                     r(i0, j0, k0, 0, 0, 0)        = d_l[fac.get_index(k0,j0,i0)][fac.num_nodes()];
                     r(i0, j0, k0, 0, height+1, 0) = d_l[fac.get_index(k0,j0,i0)][fac.num_nodes() + 1];
                } else {
                    r(i0, j0, k0, 0, 0, 0)        = d_l[fac.get_index(k0,j0,i0)][fac.num_nodes()];
                    r(i0, j0, k0, 0, 0, width+1)  = d_l[fac.get_index(k0,j0,i0)][fac.num_nodes() + 1];
                }
            }
        }
    }
    
    return array;
}

namespace py = pybind11;

PYBIND11_MODULE(PerfectMatching, m) {
    py::class_<PerfectMatching>(m, "PerfectMatching")
        .def(py::init<int, int>())
        .def("AddEdge", &PerfectMatching::AddEdge, "Add an edge to the Graph")
        .def("Solve", &PerfectMatching::Solve, "Computes a perfect matching of minimum cost.")
        .def("GetSolution", &PerfectMatching::GetSolution, "returns 1 if e is in the matching, 0 otherwise")
        .def("GetMatch", &PerfectMatching::GetMatch, "Returns the vertex that is matched with the selected vertex");

    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: cmake_example

        .. autosummary::
           :toctree: _generate

           add
           subtract
    )pbdoc";

    m.def("get_corrections", &get_corrections, R"pbdoc(
        get the corrections corresponding to the 
    )pbdoc");
    
    m.def("get_distances", &get_distances, R"pbdoc(
       get the minimum distance on the correction matrix, using BGL dijkstra.
    )pbdoc");

#ifdef VERSION_INFO
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
    m.attr("__version__") = (TOSTRING(VERSION_INFO));
#else
    m.attr("__version__") = "dev";
#endif
}
