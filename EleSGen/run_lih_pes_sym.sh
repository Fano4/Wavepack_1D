#!/bin/bash

#INCLUDING FUNCTIONS FILES
script_loc=/data1/home/stephan/Wavepack_1D/EleSGen
. ${script_loc}/wmolpro_input.sh
. ${script_loc}/extract_data.sh
. ${script_loc}/overlap_routines.sh

#DEFINING FILES LOCATION
tmp_loc=/data2/stephan/tmp_
output_loc=/data1/home/stephan/LiH_gridtest
#photoion_comp_loc=/data1/home/stephan/photoionization_coupling_comp/version_1.0
#photoion_comp_template_1=photoion_comp_template.cpp
#photoion_comp_template_2=global_vars_ref.hpp
molpro_input=LiH_
molpro_nac_input=LiH_nac_
nac_tau_output=LiH_TNAC_
nac_g_output=LiH_GNAC_
molpro_prefile=LiH_prefile
molpro_wfu_file=lih_neut_
pes_neut=LiH_neut_
pes_cat=LiH_cat_
molpro_wfu_prefile=lih_cat_
orbital_name=lih_neut
molpro_dir=/data1/apps/molpro_2012_1/bin/molpro
home_dir=/data1/home/stephan
phase_output=phase_
overlap_code_dir=/data1/home/stephan/Wavepack_1D/EleSGen/Overlap

#DEFINING PES VARIABLES

Rmin=0.8 #!!!! THESE VALUES ARE IN ANGSTROM AND NOT IN ATOMIC UNITS!!!
Rmax=21.6
grid_size=512
mH=1.007825
mLi=6.015122795
mu=$(awk "BEGIN {print 1836 * ${mH} * ${mLi} / ( ${mH} + ${mLi} ) }")

#DEFINING COMPUTATION VARIABLES
n_sym=4
n_elec=4
declare n_states_neut_sym=(8 3 3 1)
declare n_states_cat_sym=(2 1 1 0)
n_states_dicat_sym1=0
weight_neut_states_sym=("26,1,1,1,1,1,1,1" "1,1,1" "1,1,1" "1")
basis_set="{
!
! HYDROGEN       (6s,3p,1d) -> [4s,3p,1d]
! HYDROGEN       (5s)->[3s]
! HYDROGEN       (3p,1d)
! HYDROGEN       (1s)
s, H , 33.86500, 5.094790, 1.158790, 0.325840, 0.102741, 0.036000,0.01095,0.0027375
c, 1.3, 0.0254938, 0.190373, 0.852161
c, 4.4, 1
c, 5.5, 1
c, 6.6, 1
c, 7.7, 1
c, 8.8, 1
p, H , 3.0000, 0.7500, 0.1875
c, 1.1, 1
c, 2.2, 1
c, 3.3, 1
d, H , 1.0000
c, 1.1, 1
! LITHIUM       (12s,6p,3d,1f) -> [5s,4p,3d,1f]
! LITHIUM       (11s,5p)->[4s,3p]
! LITHIUM       (3d,1f)
! LITHIUM       (1s,1p)
s, LI , 900.4600, 134.4330, 30.43650, 8.626390, 2.483320, 0.303179, 4.868900, 0.856924, 0.243227, 0.0635070, 0.0243683, 0.007400
c, 1.6, 0.00228704, 0.0176350, 0.0873434, 0.2809770, 0.6587410, 0.118712
c, 7.9, 0.0933293, 0.9430450, -0.00279827
c, 10.10, 1.000000
c, 11.11, 1.000000
c, 12.12, 1.000000
p, LI , 4.868900, 0.856924, 0.243227, 0.0635070, 0.0243683, 0.007400
c, 1.3, 0.0327661, 0.1597920, 0.8856670
c, 4.4, 1.000000
c, 5.5, 1.000000
c, 6.6, 1.000000
d, LI , 0.800, 0.200, 0.050
c, 1.1, 1
c, 2.2, 1
c, 3.3, 1
f, LI , 0.150
c, 1.1, 1
}"
#6-311G++\(2df,2p\)
declare n_occ_sym=(11 4 4 1)

###############################
if [[ 0 -eq 0 ]] 
then
##############################
#FOR EACH POINT ON THE PES
echo "running PES computation"

let i=0
while [[ ${i} -lt ${grid_size} ]]
do

   echo "loop $i"
   #DEFINE THE BOND LENGTH

   R=$(awk " BEGIN { print $Rmin + $i * ( $Rmax - $Rmin ) / ($grid_size) } " )
   RH=$(awk " BEGIN { print - ( $mLi / ( $mH + $mLi ) ) * $R } " )
   RLi=$(awk " BEGIN { print ( $mH / ( $mH + $mLi ) ) * $R } " )
   dR=$(awk " BEGIN { print ( $Rmax - $Rmin ) / ($grid_size) } " )

   #CREATE A TEMP DIRECTORY
   let j=0

   fol=${tmp_loc}${j}

   while [[ -d ${fol} ]]
   do
      echo "trying temp file $j"
      let j=$j+1   
      fol=${tmp_loc}${j}
   done

   mkdir ${fol}

   echo "creating folder ${fol}"

   #CREATE MOLPRO INPUT FILES
   echo "creating molpro input file ${fol}/${molpro_prefile}.com"

#   prefiles_input #write molpro prefile input

   echo "creating molpro input file ${fol}/${molpro_input}${R}.com"

#   main_input #write main molpro input

#   echo "creating molpro input file ${fol}/${molpro_nac_input}${R}.com"

   cp ${output_loc}/${molpro_wfu_file}${R}.wfu ${fol}/
   
   NAC_input #write nac molpro input


cd ${fol}
echo "running molpro "
#${molpro_dir} -n4 -m4000M  -d ${fol} -s ${fol}/${molpro_prefile}.com
#${molpro_dir} -n4 -m4000M  -d ${fol} -s ${fol}/${molpro_input}${R}.com
${molpro_dir} -n4 -m4000M  -d ${fol} -s ${fol}/${molpro_nac_input}${R}.com
cd ${home}

#CONFIGURE PHOTOIONIZATION COMPUTATION PROGRAM

#sed "s@MO_CUBE_LOC@\"${fol}/lih_neut_orbital_\"@" ${photoion_comp_loc}/${photoion_comp_template_1} > ${photoion_comp_loc}/temp1.cpp
#sed "s@LiH_PICE_@${fol}/LiH_PICE_${R}_@" ${photoion_comp_loc}/temp1.cpp > ${photoion_comp_loc}/temp2.cpp
#sed "s@DYSON_CUBE_LOC@\"${fol}/lih_dyson_orbital_\"@" ${photoion_comp_loc}/temp2.cpp > ${photoion_comp_loc}/photoion_comp.cpp
#sed "s@molpro_output.out@${fol}/${molpro_input}${R}.out@" ${photoion_comp_loc}/${photoion_comp_template_2} > ${photoion_comp_loc}/global_vars.hpp

#RUN PHOTOIONIZATION COMPUTATION PROGRAM
#cd ${fol}
#${photoion_comp_loc}/compile.sh
#${photoion_comp_loc}/photoionization_comp.exe >> error_photoionization.txt
#cd ${home} 


#print_coordinate

#print_pes_neutral

#print_pes_cation

#print_dipole_neutral

print_nac_neutral


mv ${fol}/* ${output_loc}
rm -r ${fol}

let i=$i+1
done
############
fi
############
#RUN OVERLAP PROGRAM AND COMPUTE THE UNSCALED NON-ADIABATIC COUPLING USING FINITE DIFFERENCE METHOD.

let i=0
while [[ ${i} -lt ${grid_size} ]]
do
   let j=0

   fol=${tmp_loc}${j}

   while [[ -d ${fol} ]]
   do
      echo "trying temp file $j"
      let j=$j+1   
      fol=${tmp_loc}${j}
   done

   mkdir -p ${fol}

   echo "creating folder ${fol}"

#get_phase

#RUN OVERLAP PHASE PROGRAM AND PRINT THE PHASE VECTOR FOR THE CURRENT GEOMETRY WITH RESPECT TO THE FIRST ONE

#nac_findiff_gen
rm -r ${fol}

let i=$i+1
done
#AT LAST, MULTIPLY EACH TRANSITION QUANTITY BY ITS PHASE WITH RESPECT TO THE FIRST ELECTRONIC WAVE FUNCTION.

exit

