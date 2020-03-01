#ifndef COMPUTATION_H
#define COMPUTATION_H

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <mkl.h>
#include <omp.h>


bool Runge_kutta(wavefunction *Psi0,hamilton_matrix *H,int time_index);
bool Runge_kutta_notdH(wavefunction *Psi0,hamilton_matrix *H, int time_index);
bool t_deriv(wavefunction *Psi,hamilton_matrix *H,wavefunction *dPsi,double time_index);
bool t_deriv_matrix(hamilton_matrix *H,wavefunction **dPsi_mat,double time_index);
void propagate(wavefunction *Psi, hamilton_matrix *H,int* time_index,int num_of_loop);
bool adam_bashforth_moulton(wavefunction *dPsim4,wavefunction *dPsim3, wavefunction *dPsim2,wavefunction *dPsim1,hamilton_matrix *H,wavefunction *Psim1,int time_index);

#endif // COMPUTATION_H

