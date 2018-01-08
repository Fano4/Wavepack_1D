#!/bin/bash

#DEFINING FILES LOCATION
tmp_loc=/data2/stephan/tmp_
output_loc=/data1/home/stephan/LiH_molpro_output
photoion_comp_loc=/data1/home/stephan/photoionization_coupling_comp
photoion_comp_template_1=photoion_comp_template.cpp
photoion_comp_template_2=global_vars_ref.hpp
molpro_input=LiH_
molpro_nac_input=LiH_nac_
molpro_wfu_file=Lih_neut.wfu
molpro_dir=/data1/apps/molpro_2015_1/bin/molpro

#DEFINING PES VARIABLES

Rmin=1.63
Rmax=21.6
grid_size=4
mH=1.007825
mLi=6.015122795
dR=0.02

#DEFINING COMPUTATION VARIABLES
n_states_neut=14
n_states_cat=5
n_states_dicat=1
basis_set=6-311G++\(2df,2p\)
n_occ=18

#FOR EACH POINT ON THE PES
echo "running PES computation"

let i=0
while [[ ${i} -lt ${grid_size} ]]
do

   echo "loop $i"
   #DEFINE THE BOND LENGTH
   R=$(awk " BEGIN { print $Rmin + $i * ( $Rmax - $Rmin ) / ($grid_size-1) } " )
   RH=$(awk " BEGIN { print - ( $mLi / ( $mH + $mLi ) ) * $R } " )
   RLi=$(awk " BEGIN { print ( $mH / ( $mH + $mLi ) ) * $R } " )

   RpdR=$(awk "BEGIN{print $R + $dR }")
   RmdR=$(awk "BEGIN{print $R - $dR }")

   RpdRH=$(awk " BEGIN { print - ( $mLi / ( $mH + $mLi ) ) * $RpdR } " )
   RpdRLi=$(awk " BEGIN { print ( $mH / ( $mH + $mLi ) ) * $RpdR } " )
   RmdRH=$(awk " BEGIN { print - ( $mLi / ( $mH + $mLi ) ) * $RmdR } " )
   RmdRLi=$(awk " BEGIN { print ( $mH / ( $mH + $mLi ) ) * $RmdR } " )

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
   echo "

***,lithium_hydride
 
gprint,orbitals,civector
symmetry,nosym 
geometry={
 2
 LiH
  H          0.0000000000        0.0000000000      ${RH}
  Li         0.0000000000        0.0000000000      ${RLi}
 }
 basis=${basis_set}

{matrop;
 load,s
 print s
 }

file,2,${molpro_wfu_file},new
{casscf
maxiter, 40
closed, 0
occ, ${n_occ}
!pspace,5
config
wf,4,1,0;state, ${n_states_neut} ;
}

{cube,LiH_neut.cub,-1,125,125,125;
 orbital,occ
 BRAGG,10
 }

!{hf
!wf,3,1,1,1
!}
{casscf
maxiter,40
closed, 0
occ, ${n_occ}
!pspace,5
config
wf,3,1,1,1;state, ${n_states_cat} ;
}
!{hf
!wf,2,1,2,2
!}
{casscf
maxiter,40
closed, 0
occ, ${n_occ}
!pspace,5
config
wf,2,1,2,2;state, ${n_states_dicat} ;
}
!{hf
!wf,2,1,0,2
!}
{casscf
maxiter,40
closed, 0
occ, ${n_occ}
!pspace,5
config
wf,2,1,0,2;state, ${n_states_dicat} ;
}
---
   " > ${fol}/${molpro_input}${R}.com

   echo "creating molpro input file ${fol}/${molpro_input}${R}.com"

   echo "

***, Lithium hydride NAC

symmetry,nosym 
RH=${RH}
RLi=${RLi}
geometry={
 2
 LiH
  H          0.0000000000        0.0000000000      RH
  Li         0.0000000000        0.0000000000      RLi
 }
 basis=${basis_set}

file,2,Lih_neut.wfu
{casscf
maxiter, 40
closed, 0
occ, ${n_occ}
!pspace,5
config
wf,4,1,0;state, ${n_states_neut} ;
orbital,2140.2
}
{ci;state,${n_states_neut};noexc;
   save,6000.2;
   dm,8000.2}

RH=${RpdRH}
RLi=${RpdRLi}
{casscf
maxiter, 40
closed, 0
occ, ${n_occ}
!pspace,5
config
wf,4,1,0;state, ${n_states_neut} ;
orbital,2141.2
}
{ci;state,${n_states_neut};noexc;
   save,6001.2}
{ci;trans,6000.2,6001.2;
dm,8100.2}

RH=${RmdRH}
RLi=${RmdRLi}
{casscf
maxiter, 40
closed, 0
occ, ${n_occ}
!pspace,5
config
wf,4,1,0;state, ${n_states_neut} ;
orbital,2142.2
}
{ci;state,${n_states_neut};noexc;
   save,6002.2}
{ci;trans,6000.2,6002.2;
dm,8200.2}

   " > ${fol}/${molpro_nac_input}${R}.com

   echo "creating molpro input file ${fol}/${molpro_nac_input}${R}.com"

let m=1
while [[ $m -le $n_states_neut ]]
do
   let n=1
   while [[ $n -lt $m ]]
   do
   echo "
{ddr,2*${dR}
orbital,2140.2,2141.2,2142.2;
density,8000.2,8100.2,8200.2
STATE,${m}.1,${n}.1}
   " >> ${fol}/${molpro_nac_input}${R}.com

   let n=$n+1
done
let m=$m+1
done

cd ${fol}
echo "running molpro "
${molpro_dir} -n4 -m4000M  -d${fol} -s ${fol}/${molpro_input}${R}.com
${molpro_dir} -n4 -m4000M  -d${fol} -s ${fol}/${molpro_nac_input}${R}.com
cd

#CONFIGURE PHOTOIONIZATION COMPUTATION PROGRAM

#sed "s@MO_CUBE_LOC@\"${fol}/lih_neut_orbital_\"@" ${photoion_comp_loc}/${photoion_comp_template_1} > ${photoion_comp_loc}/temp1.cpp
#sed "s@LiH_PICE_@${fol}/LiH_PICE_${R}_@" ${photoion_comp_loc}/temp1.cpp > ${photoion_comp_loc}/temp2.cpp
#sed "s@DYSON_CUBE_LOC@\"${fol}/lih_dyson_orbital_\"@" ${photoion_comp_loc}/temp2.cpp > ${photoion_comp_loc}/photoion_comp.cpp
#sed "s@molpro_output.out@${fol}/${molpro_input}${R}.out@" ${photoion_comp_loc}/${photoion_comp_template_2} > ${photoion_comp_loc}/global_vars.hpp

#RUN PHOTOIONIZATION COMPUTATION PROGRAM
#cd ${fol}
#${photoion_comp_loc}/compile.sh
#${photoion_comp_loc}/photoionization_comp.exe
#cd 

#WRITE THE PES IN THE OUTPUT FILE

echo "${R}">>${output_loc}/coordinates.input

#GATHER ENERGIES IN MOLPRO OUTPUT

let m=1
while [[ $m -le $n_states_neut ]]
do
   grep -m1 "!MCSCF STATE ${m}.1 Energy" ${fol}/${molpro_input}${R}.out | awk '{print $5}' >> ${output_loc}/LiH_neut_${m}.1.input
   let m=$m+1
done
let m=1
while [[ $m -le $n_states_cat ]]
do
   grep -m2 "!MCSCF STATE ${m}.1 Energy" ${fol}/${molpro_input}${R}.out | tail -n1 | awk '{print $5}' >> ${output_loc}/LiH_cat_${m}.1.input
   let m=$m+1
done

#GATHER DIPOLES IN MOLPRO OUTPUT

for LL in X Y Z
do
let m=1
while [[ $m -le $n_states_neut ]]
do
   let n=1
   while [[ $n -le $m ]]
   do
      a=$(grep -m1 "<${m}.1|DM${LL}|${n}.1>" ${fol}/${molpro_input}${R}.out | awk '{print $4}')

      if [[ -z $a ]]; then
         echo "0.0" >> ${output_loc}/LiH_DM${LL}_${m}_${n}.input
      else
         echo "$a" >> ${output_loc}/LiH_DM${LL}_${m}_${n}.input
      fi
      let n=$n+1
   done
   let m=$m+1
done
done

#GATHER NAC IN MOLPRO OUTPUT

let m=1 
let q=1
while [[ $m -le $n_states_neut ]]
do
   let n=1
   while [[ $n -lt $m ]]
   do
      grep -m${q} -a5 "Transition density (R|R-DR) from     8200.2 for states"  ${fol}/${molpro_nac_input}${R}.out |tail -n1 | awk '{print $3}' >> ${output_loc}/LiH_NAC_${m}_${n}.input
      let q=$q+1
      let n=$n+1
   done
   let m=$m+1
done

rm ${fol}/*.cube
mv ${fol}/* ${output_loc}
rm -r ${fol}

let i=$i+1
done
