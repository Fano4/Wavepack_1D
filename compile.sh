#!/bin/bash
#DEFAULT=/data1/home/stephan/Wavepack_1D
#SOURCE=/data1/home/stephan/Wavepack_1D

HF5_ROOT=data1/home/stephan/libraries/cpp/hdf5
GSL_ROOT=data1/home/stephan/gsl
PICE_ROOT=data1/home/stephan/Photoionization_SAE_PW/version_3.0/photoionization_module_1.0

SOURCE=$(pwd)
DEFAULT=${SOURCE}
#HF5_ROOT=${SOURCE}/hdf5
#GSL_ROOT=${SOURCE}/gsl

if [[ -z $1 ]]
then
   DIRECTORY=$DEFAULT
else
   DIRECTORY=$1
fi

icpp=/opt/intel/composer_xe_2015.3.187/bin/intel64/icpc
${icpp} -O3 -Wall -qopenmp -I/${GSL_ROOT}/include -I/${MKLROOT}include -I/${MKLROOT}include/fftw -I/${HF5_ROOT}/include -I/${PICE_ROOT} -L/${HF5_ROOT}/lib /${HF5_ROOT}/lib/libhdf5_hl_cpp.a /${HF5_ROOT}/lib/libhdf5_cpp.a /${HF5_ROOT}/lib/libhdf5_hl.a /${HF5_ROOT}/lib/libhdf5.a -L/${GSL_ROOT}/lib -L/${MKLROOT}lib/intel64 -L/${PICE_ROOT} -Wl,-rpath -Wl,/${HF5_ROOT}/lib -mkl -lgsl -liomp5 -lhdf5 -lhdf5_cpp -lpthread -ldl -lrt -lz -DMKL_Complex16="std::complex<double>" ${DIRECTORY}/wavepack_1d.cpp -o ${DIRECTORY}/test.exe

#~/bin/wavepack.exe
