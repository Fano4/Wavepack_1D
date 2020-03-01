

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <omp.h>
#include <hdf5.h>
#include <complex>

#include "mkl.h"
#define MKL_Complex16 std::complex<double>

#include "pice_module.cpp"
#include "objects.hpp"
#include "Computation.cpp"
#include "obj_wavefunction.cpp"
#include "obj_hamilt_mat_init.cpp"
#include "obj_hamilt_mat_comp.cpp"
#include "electric_field_func.cpp"
#include "input_reader.cpp"


