#ifndef SIM_H
#define SIM_H

#include <vector>
#include <complex>

////Complex-double
//typedef struct {
//   double x,y ;
//} ds_Complex;

using ds_Complex = std::complex<double>;
using qubit_id_t = size_t;

//ds_reg struct
typedef struct {
   std::vector<ds_Complex> state; //The state vector, as an array of complex doubles on the heap.
   std::vector<int> steps;
   // int nc, nq;
   double err, sigma; // the probability of pauli errors and the std-deviation for coherant errors.
} ds_Register;

//ds_Registor Constructor
// Inputs: nq_L    - The Number of qubits in the register
//         err_L   - the probability of a Pauli-error occuring.
//         sigma_L - the standard deviation of single-gate
// Writes log to mem_regs.txt
ds_Register ds_create_register(size_t nq_L, double err_L, double sigma_L);

std::string ds_string(ds_Register& reg);

// Initialise the simulator
//void ds_initialize_simulator();
//Clear the state and steps of input
void ds_clearreg(ds_Register& reg);
//set the state element indexed by n to x+iy. no bounds check.
void ds_set_state(ds_Register& reg, size_t n, double x, double y);
//check if reg.state[n].x > 1-tol
bool ds_query_state(const ds_Register& reg, size_t n, double tol);
//Writes the state in reg to reg_out.dat
void ds_print(const ds_Register &reg);

//Add errors on each timestep until the two qubits are at the same timestep.
void ds_update(ds_Register& reg, qubit_id_t q1, qubit_id_t q2);
//bring the errors on all qubits upto the current highest timestep.
void ds_global_update(ds_Register& reg);

// Calculates the inner product of the state vectors.
double ds_inner_product(const ds_Register& reg1,const ds_Register& reg2);

//In all the single qubit gates, time is passed into lerr, and used to increment the steps vector
// for the given qubit.

//Performs an x-rot gate on qubit-q with angle theta
void ds_xrot(ds_Register& reg, qubit_id_t q, double theta, int time);
//Performs a y-rot gate on qubit-q with angle theta
void ds_yrot(ds_Register& reg, qubit_id_t q, double theta, int time);
//Performs a z-rot gate on qubit-q with angle theta
void ds_zrot(ds_Register& reg, qubit_id_t q, double theta, int time);
// Performs an X-gate on qubit q.
void ds_X(ds_Register& reg, qubit_id_t q, int time);
//Performs a Z-gate on qubit q
void ds_Z(ds_Register& reg, qubit_id_t q, int time);
//Performs an XZ-gate on qubit q (i.e. iYgate)
void ds_XZ(ds_Register& reg, qubit_id_t q, int time);
//Performs a Hadimard gate on qubit q.
void ds_Hadamard(ds_Register& reg, qubit_id_t q, int time);

//Perform required errors on the qubit q, either with random coherant errors, or random errors.
// As written will not perform coherant errors when there is 0 chance of pauli errors.
// As written will perform any errors if time=0.
void ds_lerr(ds_Register& reg, qubit_id_t q, int time);

//In all of the two-qubit gates, time is used to increment the steps vector for the given qubits and
// It is passed into the lerr function for each qubit. ds_update is called before the gate operations.

//Performs A cnot gate between qcont and qtarg
void ds_cnot(ds_Register& reg, qubit_id_t qcont, qubit_id_t qtarg, int time);
//Performs a swap gate between q1 and q2
void ds_swap(ds_Register& reg, qubit_id_t q1, qubit_id_t q2, int time);
//Performs a CZ(theta) controlled rotation on qcnt and qtarg
void ds_cphase(ds_Register& reg, qubit_id_t qcont, qubit_id_t qtarg, double theta, int time);

// measures the nq2m qubits in the list given by *lq2m, returns the results as a bit-vector (in an int)
uint64_t ds_measure(ds_Register& reg, qubit_id_t nq2m, qubit_id_t *lq2m);
// performs a measurement with a given result and returns the probability of that result.
double ds_set_measure(ds_Register& reg, qubit_id_t nq2m, qubit_id_t *lq2m, uint64_t val);
// returns the probability of returning a parity-1 measurement on the given qubits, does not perform measurement.
double ds_set_measure_parity(ds_Register& reg, size_t nq2m, qubit_id_t *lq2m);

#endif
