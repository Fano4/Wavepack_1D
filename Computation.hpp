#ifndef COMPUTATION_H
#define COMPUTATION_H

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <mkl.h>
#include <omp.h>


bool Runge_kutta(wavefunction *Psi0,hamilton_matrix *H,int time_index);
bool t_deriv(wavefunction *Psi,hamilton_matrix *H,wavefunction *dPsi,double time_index);
void propagate(wavefunction *Psi, hamilton_matrix *H,int* time_index,int num_of_loop);

#endif // COMPUTATION_H

