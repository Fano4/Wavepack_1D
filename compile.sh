#!/bin/bash
icpp=/opt/intel/composer_xe_2015.3.187/bin/intel64/icpc
${icpp} -g -qopenmp  -I/${MKLROOT}include -I/${MKLROOT}include/fftw  -L/${MKLROOT}lib/intel64 -mkl -DMKL_Complex16="std::complex<double>" -liomp5 -lpthread -ldl /data1/home/stephan/Wavepack_1D/wavepack_1d.cpp -o ~/Wavepack_1D/wavepack_1d_ceppi.exe
#~/bin/wavepack.exe
