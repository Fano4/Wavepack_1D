#!/bin/bash
#DEFAULT=/data1/home/stephan/Wavepack_1D
#SOURCE=/data1/home/stephan/Wavepack_1D

#HF5_ROOT=data1/home/stephan/libraries/cpp/hdf5
#GSL_ROOT=data1/home/stephan/gsl
#PICE_ROOT=data1/home/stephan/Photoionization_SAE_PW/version_3.0/photoionization_module_1.0

module load intel/2016b
module load HDF5/1.8.17-intel-2016b
module load GSL/2.2.1-intel-2016b

SOURCE=$(pwd)
DEFAULT=${SOURCE}
PWAPIC_ROOT=home/ulg/cpt/svdwild/Photoionization_SAE_PW/version_3.0
#HF5_ROOT=${SOURCE}/hdf5
#GSL_ROOT=${SOURCE}/gsl

if [[ -z $1 ]]
then
   DIRECTORY=$DEFAULT
else
   DIRECTORY=$1
fi

#icpp=/opt/intel/composer_xe_2015.3.187/bin/intel64/icpc
icpp=icpc
#${icpp} -I/apps/brussel/interlagos/software/zlib/1.2.11-GCCcore-6.4.0/include -I/apps/brussel/interlagos/software/Szip/2.1.1-GCCcore-6.4.0/include -I/apps/brussel/interlagos/software/imkl/2018.1.163-iimpi-2018a/mkl/include -I/apps/brussel/interlagos/software/zlib/1.2.11-GCCcore-6.4.0/include -I/apps/brussel/interlagos/software/Szip/2.1.1-GCCcore-6.4.0/include -I/apps/brussel/interlagos/software/libxml2/2.9.7-GCCcore-6.4.0/include  -I/${PWAPIC_ROOT}/photoionization_module_1.0 -O0 -Wall -mavx -ftz -fp-speculation=safe -fp-model source -fPIC -qopenmp -L/apps/brussel/interlagos/software/HDF5/1.10.1-intel-2018a/lib /apps/brussel/interlagos/software/HDF5/1.10.1-intel-2018a/lib/libhdf5_hl_cpp.a /apps/brussel/interlagos/software/HDF5/1.10.1-intel-2018a/lib/libhdf5_cpp.a /apps/brussel/interlagos/software/HDF5/1.10.1-intel-2018a/lib/libhdf5_hl.a /apps/brussel/interlagos/software/HDF5/1.10.1-intel-2018a/lib/libhdf5.a -L/apps/brussel/interlagos/software/zlib/1.2.11-GCCcore-6.4.0/lib -L/apps/brussel/interlagos/software/Szip/2.1.1-GCCcore-6.4.0/lib -L/apps/brussel/interlagos/software/icc/2018.1.163-GCC-6.4.0-2.28/lib/intel64 -L/apps/brussel/interlagos/software/imkl/2018.1.163-iimpi-2018a/lib -L/apps/brussel/interlagos/software/imkl/2018.1.163-iimpi-2018a/mkl/lib/intel64 -L/apps/brussel/interlagos/software/imkl/2018.1.163-iimpi-2018a/lib -L/apps/brussel/interlagos/software/zlib/1.2.11-GCCcore-6.4.0/lib -L/apps/brussel/interlagos/software/Szip/2.1.1-GCCcore-6.4.0/lib -L/apps/brussel/interlagos/software/libxml2/2.9.7-GCCcore-6.4.0/lib -L/${PWAPIC_ROOT}/photoionization_module_1.0 -mkl -lgsl -liomp5 -lhdf5 -lhdf5_cpp -lpthread -ldl -lrt -lsz -lz -DMKL_Complex16="std::complex<double>" -Wl,-rpath -Wl,/apps/brussel/interlagos/software/HDF5/1.10.1-intel-2018a/lib ${DIRECTORY}/wavepack_1d.cpp -o ${DIRECTORY}/wavepack_1d.exe
${icpp} -I/apps/brussel/interlagos/software/zlib/1.2.8-intel-2016b/include -I/apps/brussel/interlagos/software/Szip/2.1-intel-2016b/include -I/apps/brussel/ivybridge/software/imkl/11.3.3.210-iimpi-2016b/mkl/include -I/apps/brussel/interlagos/software/zlib/1.2.8-intel-2016b/include -I/apps/brussel/interlagos/software/Szip/2.1-intel-2016b/include -I/${PWAPIC_ROOT}/photoionization_module_1.0 -fPIC -O2 -qopenmp -xHost -ftz -fp-speculation=safe -fp-model source -L/apps/brussel/interlagos/software/HDF5/1.8.17-intel-2016b/lib /apps/brussel/interlagos/software/HDF5/1.8.17-intel-2016b/lib/libhdf5_hl_cpp.a /apps/brussel/interlagos/software/HDF5/1.8.17-intel-2016b/lib/libhdf5_cpp.a /apps/brussel/interlagos/software/HDF5/1.8.17-intel-2016b/lib/libhdf5_hl.a /apps/brussel/interlagos/software/HDF5/1.8.17-intel-2016b/lib/libhdf5.a -L/apps/brussel/interlagos/software/zlib/1.2.8-intel-2016b/lib -L/apps/brussel/interlagos/software/Szip/2.1-intel-2016b/lib -L/apps/brussel/ivybridge/software/icc/2016.3.210-GCC-5.4.0-2.26/lib/intel64 -L/apps/brussel/ivybridge/software/imkl/11.3.3.210-iimpi-2016b/lib -L/apps/brussel/ivybridge/software/imkl/11.3.3.210-iimpi-2016b/mkl/lib/intel64 -L/apps/brussel/ivybridge/software/imkl/11.3.3.210-iimpi-2016b/lib -L/apps/brussel/interlagos/software/zlib/1.2.8-intel-2016b/lib -L/apps/brussel/interlagos/software/Szip/2.1-intel-2016b/lib -L/${PWAPIC_ROOT}/photoionization_module_1.0 -mkl -lgsl -liomp5 -lhdf5 -lhdf5_cpp -lrt -lsz -lz -ldl -lpthread -DMKL_Complex16="std::complex<double>" -Wl,-rpath -Wl,/apps/brussel/interlagos/software/HDF5/1.8.17-intel-2016b/lib ${DIRECTORY}/wavepack_1d.cpp -o ${DIRECTORY}/wavepack_1d.exe
#icpp=/opt/intel/composer_xe_2015.3.187/bin/intel64/icpc
#${icpp} -O3 -g -Wall -qopenmp -I/${GSL_ROOT}/include -I/${MKLROOT}include -I/${MKLROOT}include/fftw -I/${HF5_ROOT}/include -I/${PICE_ROOT} -L/${HF5_ROOT}/lib /${HF5_ROOT}/lib/libhdf5_hl_cpp.a /${HF5_ROOT}/lib/libhdf5_cpp.a /${HF5_ROOT}/lib/libhdf5_hl.a /${HF5_ROOT}/lib/libhdf5.a -L/${GSL_ROOT}/lib -L/${MKLROOT}lib/intel64 -L/${PICE_ROOT} -Wl,-rpath -Wl,/${HF5_ROOT}/lib -mkl -lgsl -liomp5 -lhdf5 -lhdf5_cpp -lpthread -ldl -lrt -lz -DMKL_Complex16="std::complex<double>" ${DIRECTORY}/wavepack_1d.cpp -o ${DIRECTORY}/test.exe

#~/bin/wavepack.exe
