#!/bin/bash

#The first argument is the pump probe delay in femtosecond

if [[ -z $1 ]]
then
   echo "Error. pump-probe delay not parametrized. exit"
   exit
fi
#echo $(awk "BEGIN {print $1/0.02418884}")
if [[ -z $2 ]]
then
   echo "Error. INPUT DIR NOT SPECIFIED. exit"
   exit
fi
if [[ -z $3 ]]
then
   echo "Error. OUTPUT DIR NOT SPECIFIED. exit"
   exit
fi

INPUT_DIR=$2
OUTPUT_DIR=$3
#OUTPUT_DIR=/home/ulg/cpt/svdwild/LiH_wvpck_output

#Input files rootnames
neut_potential_root=${INPUT_DIR}/Int_pes_
cat_potential_root=${INPUT_DIR}/Int_pes_cat_
neut_dipole_root=${INPUT_DIR}/LiH_neut_
cat_dipole_root=${INPUT_DIR}/LiH_cat_
neut_nac_root=${INPUT_DIR}/LiH_NAC_
neut_pice_root=${INPUT_DIR}/LiH_pice_data_newdyson.h5

#Output files rootnames
output_root=${OUTPUT_DIR}/output_${1}.log
continuum_root=${OUTPUT_DIR}/continuum_${1}_
neut_wfu_root=${OUTPUT_DIR}/neut_wf_${1}_
spectrum_root=${OUTPUT_DIR}/spectrum_${1}.txt
mfpad_root=${OUTPUT_DIR}/mfpad_${1}.txt
cs_root=${OUTPUT_DIR}/cs_${1}_
dist_root=${OUTPUT_DIR}/sphere_dist.txt

#Simulation parameters

grid_size=512
small_grid_size=128
n_states_neut=10
n_states_cat=1
n_points_sphere=128
n_k=56
k_sample=35
kmin=0.0001
kmax=1.5
xmin=$(awk 'BEGIN {print 0.8/0.529}')
xmax=$(awk 'BEGIN {print 21.6/0.529}')
mass=$(awk 'BEGIN {print (1.007825*7.01600455/(1.007825+7.01600455))*1836}')
total_time=$(awk 'BEGIN {print 230/0.02418884}')
h=$(awk 'BEGIN {print 0.0004/0.02418884}')
efield_thresh=0.0001
pot_vec_thresh=0.005

#Pump and probe parameters

pump_strength=0.01
pump_origin=$(awk 'BEGIN {print 5.0/0.02418884}')
pump_sigma=$(awk 'BEGIN {print 0.5945252602016133/0.02418884}') #$(awk 'BEGIN {print 1.1890505204/0.02418884}')
pump_energy=$(awk 'BEGIN {print 6.01136/27.211}')
pump_CEP=0.0
pprobe_delay=$(awk "BEGIN {print $1/0.02418884}")
probe_strength=0.00 #0.005
probe_sigma=$(awk 'BEGIN {print 0.30800/0.02418884}')
probe_energy=$(awk 'BEGIN {print 19.53692/27.211}')
probe_CEP=0.0


#Write the input using the parameters defined above

echo "
${neut_potential_root}
${cat_potential_root}
${neut_dipole_root}
${cat_dipole_root}
${neut_nac_root}
${neut_pice_root}
${output_root}
${continuum_root}
${neut_wfu_root}
${spectrum_root}
${mfpad_root}
${cs_root}
${dist_root}
${grid_size}
${small_grid_size}
${n_states_neut}
${n_states_cat}
${n_points_sphere}
${n_k}
${k_sample}
${kmin}
${kmax}
${xmin}
${xmax}
${mass}
${total_time}
${h}
${efield_thresh}
${pot_vec_thresh}
${pump_strength}
${pump_origin}
${pump_sigma}
${pump_energy}
${pump_CEP}
${pprobe_delay}
${probe_strength}
${probe_sigma}
${probe_energy}
${probe_CEP}
">${OUTPUT_DIR}/wvpck_input.in

