bin_PROGRAMS = wavepack_1d



wavepack_1d_SOURCES= wavepack_1d.cpp 

include_HEADERS = *.hpp

AM_CXXFLAGS= -O3 -qopenmp -v -DMKL_Complex16="std::complex<double>" -I. -I${MKLROOT}include -I/${HF5_ROOT}/include -I/${PICE_ROOT} -Wl,-rpath -Wl,/${HF5_ROOT}/lib 

AM_LDFLAGS=-L/${HF5_ROOT}/lib /${HF5_ROOT}/lib/libhdf5_hl_cpp.a /${HF5_ROOT}/lib/libhdf5_cpp.a /${HF5_ROOT}/lib/libhdf5_hl.a /${HF5_ROOT}/lib/libhdf5.a -L/${MKLROOT}lib/intel64 -L/${PICE_ROOT} -mkl -liomp5 -lhdf5 -lhdf5_cpp -lrt -ldl -lpthread 
