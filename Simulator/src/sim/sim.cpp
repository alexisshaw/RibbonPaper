#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <bit>
#include <algorithm>
#include <numeric>
#include <numbers>
#include "sim/sim.h"
#include "sim/norm.h"
#include "boost/container/small_vector.hpp"

using namespace std;

//rotate the state vector (zPtr1, zptr2)^T by theta around the x axis
void ds_esigx(ds_Complex &zPtr1, ds_Complex &zPtr2, double theta);
//rotate the state vector (zPtr1, zptr2)^T by theta around the y axis
void ds_esigy(ds_Complex &zPtr1, ds_Complex &zPtr2, double theta);
//rotate the state vector (zPtr1, zptr2)^T by theta around the z axis
void ds_esigz(ds_Complex &zPtr1, ds_Complex &zPtr2, double theta);
//rotate the state vector (zPtr1, zptr2)^T by a unitary defined by rotations, alpha,beta,theta.
void ds_unitary(ds_Complex &zPtr1, ds_Complex &zPtr2,
                double alpha, double beta, double theta);

ds_Register ds_create_register(size_t nq_L, double err_L, double sigma_L) {
    ds_Register reg;

    reg.err = err_L;
    reg.sigma = sigma_L;
    // reg.nc = 1 << nq_L;
    // reg.nq = nq_L;

    reg.state = std::vector<ds_Complex>(1 << nq_L, {0.0,0.0});
    reg.steps = std::vector<int>(nq_L, 0);

    return reg;
}

void ds_clearreg(ds_Register& reg){
    std::fill(begin(reg.state), end(reg.state), 0.0);
    std::fill(begin(reg.steps), end(reg.steps), 0);
}

void ds_set_state(ds_Register& reg, size_t n, double x, double y){
   reg.state[n] = complex(x,y);
}

bool ds_query_state(ds_Register& reg, size_t n, double tol){
    return (reg.state[n].real() >= 1-tol);
}

std::string ds_string(ds_Register& reg) {
    std::stringstream outbuf{};
    for(size_t i = 0; i < reg.state.size(); i++)
        if(std::norm(reg.state[i]) > 1E-16)
            outbuf << i << " " << reg.state[i] << '\n';
    return outbuf.str();
}

void ds_update(ds_Register& reg, qubit_id_t q1, qubit_id_t q2){
    int &sq1 = reg.steps[q1], &sq2 = reg.steps[q2];
    qubit_id_t min_index = (sq1 < sq2)? q1: q2,  max_index = (sq1 < sq2)? q2: q1;

    for (;reg.steps[min_index] < reg.steps[max_index]; reg.steps[min_index]++) {
      ds_lerr(reg, min_index, 1);
   }
}

void ds_global_update(ds_Register& reg){
   int max_steps = *(std::max_element(begin(reg.steps), end(reg.steps)));

   for (size_t q=0; q<reg.steps.size(); q++)
      for (int i=reg.steps[q]; i<max_steps; i++)
          ds_lerr(reg, q, 1);

   std::fill(begin(reg.steps), end(reg.steps), max_steps);
}

double ds_inner_product(const ds_Register& reg1, const ds_Register& reg2){
   ds_Complex temp = 0.0;

   for (size_t i=0; i<reg1.state.size(); i++)
       temp += reg1.state[i] * conj(reg2.state[i]);

   return std::norm(temp);
}

void ds_esigx(ds_Complex &zPtr1, ds_Complex &zPtr2, double theta){
   double c = cos(theta/2), s = sin(theta/2);
    ds_Complex z1 = zPtr1, z2 = zPtr2;
   zPtr1 = complex(  c*z1.real() - s*z2.imag() , c*z1.imag() + s*z2.real());
   zPtr2 = complex( -s*z1.imag() + c*z2.real() , s*z1.real() + c*z2.imag());
}

void ds_esigy(ds_Complex &zPtr1, ds_Complex &zPtr2, double theta){
    double c = cos(theta/2), s = sin(theta/2);
    ds_Complex z1 = zPtr1, z2 = zPtr2;
    zPtr1 = complex( c*z1.real() + s*z2.real() , c*z1.imag() + s*z2.imag());
    zPtr2 = complex(-s*z1.real() + c*z2.real() ,-s*z1.imag() + c*z2.imag());
}

void ds_esigz(ds_Complex &zPtr1, ds_Complex &zPtr2, double theta){
    double c = cos(theta/2), s = sin(theta/2);
    ds_Complex z1 = zPtr1, z2 = zPtr2;
    zPtr1 = complex(c*z1.real() - s*z1.imag() , c*z1.imag() + s*z1.real());
    zPtr2 = complex(c*z2.real() + s*z2.imag() , c*z2.imag() - s*z2.real());
}

void ds_unitary(ds_Complex &zPtr1, ds_Complex &zPtr2,
             double alpha, double beta, double theta){
   double c = cos(theta/2), s = sin(theta/2);
   double c11 = cos(+alpha/2+beta/2);
   double s11 = sin(+alpha/2+beta/2);
   double c12 = cos(+alpha/2-beta/2);
   double s12 = sin(+alpha/2-beta/2);
   ds_Complex z1 = zPtr1, z2 = zPtr2;

   zPtr1 = complex(c*(c11*z1.real()-s11*z1.imag()) + s*(c12*z2.real()-s12*z2.imag()),
                    c*(c11*z1.imag()+s11*z1.real()) + s*(c12*z2.imag()+s12*z2.real()));
   zPtr2 = complex(-s*(c12*z1.real()+s12*z1.imag()) + c*(c11*z2.real()+s11*z2.imag()),
                      -s*(c12*z1.imag()-s12*z1.real()) + c*(c11*z2.imag()-s11*z2.real()));
}

/* Function used when performing a single qubit operation. n runs from 0
   to nc/2 - 1 and simply keeps track of how far through the array the
   calculation is. q runs from 0 to nq - 1 and denotes which qubit the
   operation is being performed on. The function outputs the next pair of
   indices to apply the desired transformation to.

     returns: nth pair of indices *iPtr and *jPtr, where in state *iPtr
              qubit q is not set and in state *jPtr qubit q is set
     e.g. for 3 qubits: n=0, q=1 -> *iPtr=0 (000) *jPtr=1 (001)
     e.g. for 3 qubits: n=0, q=2 -> *iPtr=0 (000) *jPtr=2 (010)
     e.g. for 3 qubits: n=1, q=3 -> *iPtr=1 (001) *jPtr=1 (101)
*/
std::pair<uint64_t, uint64_t> ds_one_qubit_indices(uint64_t n, qubit_id_t q){
   /* l is used as a mask to calculate nl --- the low part of n */
   /* calculate mask of 1s for section of n less significant than bit q */
   uint64_t l = (1 << q) - 1;

   /* calculate section of n less significant than bit q */
   uint64_t nl = n & l;

   /* remove section of n less significant than bit q and shift n higher for
      later insertion of bit q */
   n = (n - nl) << 1;

   return std::make_pair(n + nl, n + nl + l + 1);
}

std::tuple<uint64_t, uint64_t, uint64_t> ds_2q_indices_common(uint64_t n, qubit_id_t qm, qubit_id_t ql){
    /* calculate mask of 1s for section of n less significant than bit ql */
    uint64_t l = (1 << ql) - 1;
    /* calculate mask of 1s for section of n less significant than bit qm */
    uint64_t m = (1 << qm) - 1;

    /* calculate section of n less significant than bit ql */
    uint64_t nl = n & l;

    /* remove section of n less significant than bit ql and shift n higher for
       later insertion of bit q2 */
    n = (n - nl) << 1;

    /* calculate section of n less significant than bit qm but more significant
       than bit ql */
    uint64_t nm = n & m;

    /* remove section of n less significant than bit qm and shift n higher for
       later insertion of bit qm. Add nm and nl to complete preparation of n */
    n = ((n - nm) << 1) + nm + nl;

    //    n = ((n << 2) & ((~m)<<1)) | ((n << 1) & (m & ((~l)<<1))) | (n & l);

    return std::make_tuple(n, m + 1, l + 1);
}

std::tuple<uint64_t, uint64_t, uint64_t, uint64_t> ds_two_qubit_indices(uint64_t n_in, qubit_id_t q1, qubit_id_t q2){
   auto [n,m,l] = ds_2q_indices_common(n_in, std::max(q1, q2), std::min(q1, q2));

   return std::make_tuple(n, n+l, n+m, n+m+l);
}

void ds_Hadamard(ds_Register& reg, qubit_id_t q, int time){
    constexpr double ds_root2_2 = 1/std::numbers::sqrt2;

    for (size_t i=0; i<reg.state.size()/2; i++) {
        auto [j,k] = ds_one_qubit_indices(i, q);

        ds_Complex z1 = reg.state[j];
        ds_Complex z2 = reg.state[k];

        reg.state[j] = ds_root2_2*(z1 + z2);
        reg.state[k] = ds_root2_2*(z1 - z2);
    }

    ds_lerr(reg, q, time);
    reg.steps[q]+=time;
}

void ds_xrot(ds_Register& reg, qubit_id_t q, double theta, int time){
   for (size_t i=0; i<reg.state.size()/2; i++) {
       auto [j,k] = ds_one_qubit_indices(i, q);
       ds_esigx(reg.state[j], reg.state[k], theta);
   }

   ds_lerr(reg, q, time);
   reg.steps[q]+=time;
}

void ds_yrot(ds_Register& reg, qubit_id_t q, double theta, int time){
   for (size_t i=0; i<reg.state.size()/2; i++) {
       auto [j,k] = ds_one_qubit_indices(i, q);
       ds_esigy(reg.state[j], reg.state[k], theta);
   }

   ds_lerr(reg, q, time);
   reg.steps[q]+=time;
}

void ds_zrot(ds_Register& reg, qubit_id_t q, double theta, int time){
   for (size_t i=0; i<reg.state.size()/2; i++) {
       auto [j,k] = ds_one_qubit_indices(i, q);
      ds_esigz(reg.state[j], reg.state[k], theta);
   }

   ds_lerr(reg, q, time);
   reg.steps[q]+=time;
}

void ds_X(ds_Register& reg, qubit_id_t q, int time){
   for (size_t i=0; i<reg.state.size()/2; i++) {
       auto [j,k] = ds_one_qubit_indices(i, q);
       std::swap(reg.state[j], reg.state[k]);
   }

   ds_lerr(reg, q, time);
   reg.steps[q]+=time;
}

void ds_Z(ds_Register& reg, qubit_id_t q, int time){
   for (size_t i=0; i<reg.state.size()/2; i++) {
       auto [j,k] = ds_one_qubit_indices(i, q);
      reg.state[k] = -reg.state[k];
   }

   ds_lerr(reg, q, time);
   reg.steps[q]+=time;
}

void ds_XZ(ds_Register& reg, qubit_id_t q, int time){
   for (size_t i=0; i<reg.state.size()/2; i++) {
       auto [j,k] = ds_one_qubit_indices(i, q);
      ds_Complex z = reg.state[j];
      reg.state[j] = -reg.state[k];
      reg.state[k] = z;
   }

   ds_lerr(reg, q, time);
   reg.steps[q]+=time;
}

void ds_lerr(ds_Register &reg, qubit_id_t q, int time){
    if (time == 0)
        return;
    if ((reg.err < 1E-14) && (reg.sigma <  1E-14))
        return;

    if (reg.err > 0){
        double p = ds_uniform();
        if (p < reg.err / 3){
            ds_X(reg, q, 0);
            // std::cout << "x" << q << ' ';
        } else if (p < 2 * reg.err / 3){
            ds_Z(reg, q, 0);
            // std::cout << "z" << q << ' ';
        } else if (p < reg.err){
            ds_XZ(reg, q, 0);
            // std::cout << "y" << q << ' ';
        }
    }

    if (reg.sigma > 1E-14){
        double alpha = reg.sigma * ds_norm();
        double beta = reg.sigma * ds_norm();
        double theta = reg.sigma * ds_norm();
        for (size_t i = 0; i < reg.state.size() / 2; i++){
            auto [j,k] = ds_one_qubit_indices(i, q);
            ds_unitary(reg.state[j], reg.state[k], alpha, beta, theta);
        }
    }
}

/* See comments for ds_one_qubit_indices(). */
std::pair<uint64_t, uint64_t> ds_controlled_indices(uint64_t n_in, qubit_id_t qcont, qubit_id_t qtarg){
    qubit_id_t qm = std::max(qcont, qtarg), ql = std::min(qcont, qtarg);

    auto [n,m,l] = ds_2q_indices_common(n_in, qm, ql);

    return std::make_pair((qcont > qtarg )? n + m : n+l, n+m+l);
}

/* Kane compatible if |qcont-qtarg| = 1 */
void ds_cnot(ds_Register& reg, qubit_id_t qcont, qubit_id_t qtarg, int time){
   ds_update(reg, qcont, qtarg);

   for (size_t i=0; i<reg.state.size()/4; i++) {
       auto [j,k] = ds_controlled_indices(i, qcont, qtarg);
      /* Note that the ds_toffoli gate does not work with pseudo ds_cnot gates
         ie : you cant use ds_esigy(reg.state+j, reg.state+k, Pi) */
      std::swap(reg.state[j], reg.state[k]);
   }

   ds_lerr(reg, qcont, time);
   reg.steps[qcont]+=time;
   ds_lerr(reg, qtarg, time);
   reg.steps[qtarg]+=time;
}


/* Kane compatible if |q-x| = 1 */

void ds_cphase(ds_Register& reg, qubit_id_t qcont, qubit_id_t qtarg, double theta, int time){
   ds_update(reg, qcont, qtarg);

   for (size_t i=0; i<reg.state.size()/4; i++) {
       auto [j,k] = ds_controlled_indices(i, qcont, qtarg);
      reg.state[k] = std::polar(1.0, theta) * reg.state[k];
   }

   ds_lerr(reg, qcont, time);
   reg.steps[qcont]+=time;
   ds_lerr(reg, qtarg, time);
   reg.steps[qtarg]+=time;
}

/* See comments for ds_one_qubit_indices(). */
std::pair<uint64_t, uint64_t> ds_swap_indices(uint64_t n_in, qubit_id_t q1, qubit_id_t q2){
    auto [n,m,l] = ds_2q_indices_common(n_in, std::max(q1, q2), std::min(q1, q2));

    return std::make_pair(n+l, n+m);
}

/* Kane compatible if |q1-q2| = 1 */
void ds_swap(ds_Register& reg, qubit_id_t q1, qubit_id_t q2, int time){
    ds_update(reg, q1, q2);

    for (size_t i=0; i<reg.state.size()/4; i++) {
        auto [j,k] = ds_swap_indices(i, q1, q2);
        std::swap(reg.state[j],reg.state[k]);
    }

    ds_lerr(reg, q1, time);
    reg.steps[q1]+=time;
    ds_lerr(reg, q2, time);
    reg.steps[q2]+=time;
}


//Measurement helper function to get an array of masks, for a given list of qubits to measure;
inline boost::container::small_vector<uint64_t, 8> get_measurement_masks(size_t nq2m, const qubit_id_t* lq2m){
    boost::container::small_vector<uint64_t, 8> masks(nq2m, 1);
    for (size_t i=0; i<nq2m; i++) {
        masks[i] = 1;
        masks[i] <<= lq2m[i];
    }
    return masks;
}

//Given a set of measurement masks representing qubits, return an array of probabilities for each possible measurement outcome.
inline boost::container::small_vector<double, 256> get_measurement_probabilities(const ds_Register& reg,
                                                                                 const boost::container::small_vector<uint64_t, 8>& masks){
    boost::container::small_vector<double, 256> probabilities(1<<masks.size(), 0.0);
    //probabilities = (double *)calloc(1 << nq2m, sizeof(double));

    //std::cout << reg.state << std::endl;
    for (size_t i=0; i<reg.state.size(); i++) {
        /* reduce array index i */
        int j = 0;
        for (auto mask : masks) {
            j <<= 1;
            if (i & mask) {
                j++;
            }
        }
        probabilities[j] += std::norm(reg.state[i]);
    }
    return probabilities;
}

//After choosing a given measurement value, and given the norm associated with the probability of that outcome,
//Correct the state so that the given
void measurement_fixup(ds_Register& reg, size_t imeas, const boost::container::small_vector<uint64_t, 8>& masks, double norm){
    size_t mask = std::accumulate(masks.begin(), masks.end(), 0);

    size_t spread = 0;
    for (size_t i = 1, j = 0; j < masks.size(); j++) {
        if (imeas & i)
            spread += masks[masks.size() - 1 - j];
        i <<= 1;
    }

    for (size_t i=0; i<reg.state.size(); i++) {
        reg.state[i] *=  ((i & mask) == spread)? norm: 0.0;
    }
}

/* nq2m : number of qubits to ds_measure
   lq2m : list of qubits to ds_measure (sorted from largest to smallest) */
uint64_t ds_measure(ds_Register& reg, size_t nq2m, qubit_id_t *lq2m){
    ds_global_update(reg);

    /* create array of single bit masks for reducing array index */
    boost::container::small_vector<uint64_t, 8> masks = get_measurement_masks(nq2m, lq2m);

    /* create array of possible measured values */
    boost::container::small_vector<double, 256> probabilities = get_measurement_probabilities(reg, masks);


    /* random value between 0 and 1 */
    double r = ds_uniform();

    // Find the corresponding value
    size_t imeas = 0;
    for(double p = probabilities[0]; (p < r) && (imeas < probabilities.size()); p += probabilities[imeas]){
        imeas++;
    }

    double norm = 1/sqrt(probabilities[imeas]);

    measurement_fixup(reg, imeas, masks, norm);

    return int(imeas);
}

/*------------set/ds_measure----------------------------*/

/* nq2m : number of qubits to measure
   lq2m : list of qubits to measure (sorted from largest to smallest)
   val  : value you wish to observe */

double ds_set_measure(ds_Register& reg, size_t nq2m, qubit_id_t *lq2m, uint64_t val){
  /* probably not necessary but simplifies things */
  ds_global_update(reg);

  /* create array of single bit masks for reducing array index */
  boost::container::small_vector<uint64_t, 8> masks = get_measurement_masks(nq2m, lq2m);

  /* create array of possible measured values */
  boost::container::small_vector<double, 256> probabilities = get_measurement_probabilities(reg, masks);

  //Get the probability of the selected value
  double p = probabilities[val];

  double norm = (probabilities[val] == 0) ? 0 : 1 / sqrt(probabilities[val]);

  measurement_fixup(reg, val, masks, norm);

  return p;
}

double ds_set_measure_parity(ds_Register& reg, size_t nq2m, qubit_id_t *lq2m){
    /* probably not necessary but simplifies things */
    ds_global_update(reg);

    /* create array of single bit masks for reducing array index */
    boost::container::small_vector<uint64_t, 8> masks = get_measurement_masks(nq2m, lq2m);

    /* create array of possible measured values */
    boost::container::small_vector<double, 256> probabilities = get_measurement_probabilities(reg, masks);

    double parity_1_p = 0.0;
    for(size_t i = 0; i < probabilities.size(); i++){
        //std::cout << i  << ' ' << std::popcount(i) << ' ' << probabilities[i] << ' ';
        if(std::popcount(i) %2 == 1){
            parity_1_p += probabilities[i];
        }
    }
    return parity_1_p;
}
