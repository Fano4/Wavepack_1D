#!/bin/bash
icpp=/opt/intel/composer_xe_2015.3.187/bin/intel64/icpc
${icpp} -g -qopenmp  -I/${MKLROOT}include -I/${MKLROOT}include/fftw  -L/${MKLROOT}lib/intel64 -mkl -liomp5 -lpthread -ldl  /data1/home/stephan/Wavepack_1D/EleSGen/Overlap/wf_overlap.cpp -o /data1/home/stephan/Wavepack_1D/EleSGen/Overlap/wf_overlap.exe
