#!/bin/bash

#INCLUDING FUNCTIONS FILES
. /data1/home/stephan/LiH_molpro_output/run_dir/wmolpro_input.sh

#DEFINING FILES LOCATION
tmp_loc=/data2/stephan/tmp_
output_loc=/data1/home/stephan/LiH_gridtest
#photoion_comp_loc=/data1/home/stephan/photoionization_coupling_comp/version_1.0
#photoion_comp_template_1=photoion_comp_template.cpp
#photoion_comp_template_2=global_vars_ref.hpp
molpro_input=LiH_
molpro_nac_input=LiH_nac_
molpro_prefile=LiH_prefile
molpro_wfu_file=lih_neut_
molpro_wfu_prefile=lih_cat_
molpro_dir=/data1/apps/molpro_2012_1/bin/molpro
home_dir=/data1/home/stephan

#DEFINING PES VARIABLES

Rmin=0.8 #!!!! THESE VALUES ARE IN ANGSTROM AND NOT IN ATOMIC UNITS!!!
Rmax=21.6
grid_size=512
mH=1.007825
mLi=6.015122795

#DEFINING COMPUTATION VARIABLES
n_sym=4
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

   prefiles_input #write molpro prefile input

   echo "creating molpro input file ${fol}/${molpro_input}${R}.com"

   main_input #write main molpro input

   echo "creating molpro input file ${fol}/${molpro_nac_input}${R}.com"

 #  NAC_input #write nac molpro input


cd ${fol}
echo "running molpro "
${molpro_dir} -n4 -m4000M  -d ${fol} -s ${fol}/${molpro_prefile}.com
${molpro_dir} -n4 -m4000M  -d ${fol} -s ${fol}/${molpro_input}${R}.com
# ${molpro_dir} -n4 -m4000M  -d ${fol} -s ${fol}/${molpro_nac_input}${R}.com
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

#WRITE THE COORDINATE IN THE OUTPUT FILE

echo "${R}">>${output_loc}/coordinates.input

#GATHER ENERGIES IN MOLPRO OUTPUT

let k=0
let ts=0
while [[ $k -lt ${n_sym} ]]
do
   let kp=$k+1
let m=1
while [[ $m -le ${n_states_neut_sym[$k]} ]]
do
   let ts=$ts+1
   echo "seeking for energy of neutral state ${m}.${kp}"
   grep -m1 "!MCSCF STATE ${m}.${kp} Energy" ${fol}/${molpro_input}${R}.out | awk '{print $5}' >> ${output_loc}/LiH_neut_${ts}.input
   let m=$m+1
done
let k=$k+1
done

let k=0
let ts=0
while [[ ${k} -lt ${n_sym} ]]
do
let kp=$k+1
let m=1
while [[ $m -le ${n_states_cat_sym[$k]} ]]
do
   let ts=$ts+1
   echo "seeking for energy of cation state ${m}.${kp}"
   grep -m2 "!MCSCF STATE ${m}.${kp} Energy" ${fol}/${molpro_input}${R}.out | tail -n1 | awk '{print $5}' >> ${output_loc}/LiH_cat_${ts}.input
   let m=$m+1
done
let k=$k+1
done

#GATHER DIPOLES AND NAC IN MOLPRO OUTPUT

#loop between the different componenets of the vector
for LL in X Y Z
do
  let o=0 #symmetry of the initial state
  #loop over symmetries of the initial state
  while [[ $o -lt $n_sym ]]
  do
    let op=$o+1 #symmetry of initial state in molpro output

    let p=$o #symmetry of the final state
    #loop over symmetries of the final state
    while [[ $p -lt $n_sym ]]
    do
       let pp=$p+1 #symmetry of the final state in molpro output

       if [[ $o -eq $p ]] # If the symmetry of the initial state is the same as the symmetry of the final state we are in a diagonal block of the matrix. we thus consider only half the elements
        then
        let ts1=0 #total index of the initial state
        let temp1=0 
        #For a given symmetry, the total index of a state is the sum over all preivous symmetries of the number of states in each symmetry
    ##############
        while [[ $temp1 -lt ${o} ]]
        do
           let ts1=${ts1}+${n_states_neut_sym[$temp1]} 
           let temp1=${temp1}+1
        done
    ##############
        let m=1 #initial state index in symetry o
        #loop over states indexes of the initial state in symmetry o
        while [[ $m -le ${n_states_neut_sym[$o]} ]]
        do
           let n=$m #state index of the final state in symmetry o
           let ts1=${ts1}+1 

           let ts2=${ts1}-1 #total index of the final state

           #loop over the states indexes of the final state in symmetry p
           while [[ $n -le ${n_states_neut_sym[$p]} ]]
           do
             let ts2=${ts2}+1
             echo "seeking for dipole moment of neutral  ${n}.${pp}-${m}.${op} = ${ts2}-${ts1}"
             a=$(grep -m1 "<${n}.${pp}|DM${LL}|${m}.${op}>" ${fol}/${molpro_input}${R}.out | awk '{print $4}')

             #if there is no instance of the sought dipole, put it to zero
             if [[ -z $a ]]; then
               echo "0.0" >> ${output_loc}/LiH_DM${LL}_${ts2}_${ts1}.input
             else
               echo "$a" >> ${output_loc}/LiH_DM${LL}_${ts2}_${ts1}.input
             fi

             #Once for all orientations, seek for the NAC. There is no NAC between states of different symmetries
             if [[ $LL -eq X ]]; then
               a=$(grep "SA-MC NACME FOR STATES ${n}.${pp} - ${m}.${op}"  ${fol}/${molpro_nac_input}${R}.out)
               if [[ -z $a ]]; then
                  let n=$n+1
                  continue
               else
                  var1=$(grep -a4 "SA-MC NACME FOR STATES ${n}.${pp} - ${m}.${op}"  ${fol}/${molpro_nac_input}${R}.out |tail -n1|awk '{print $4}')
                  var2=$(grep -a5 "SA-MC NACME FOR STATES ${n}.${pp} - ${m}.${op}"  ${fol}/${molpro_nac_input}${R}.out |tail -n1|awk '{print $4}')
                  res=$(awk "BEGIN {print ( $mLi / ( $mH + $mH ) ) * $var2 - ( $mLi / ( $mLi + $mH ) ) * $var1 }")
                 echo "$res" >> ${output_loc}/LiH_NAC_${ts2}_${ts1}.input
               fi
             fi
             let n=$n+1
          done
          let m=$m+1
        done
     else #if initial sym is not the same as final sym, we are in an off-diagonal block of the matrix. we have ro consider all elements
        let m=1
        let ts1=0
        let temp1=0
        ##############
        while [[ $temp1 -lt ${o} ]]
        do
          let ts1=${ts1}+${n_states_neut_sym[$temp1]} 
          let temp1=${temp1}+1
        done
        ##############

        while [[ $m -le ${n_states_neut_sym[$o]} ]]
        do
          let ts1=${ts1}+1
          let n=1
          let ts2=0
          let temp1=0
          ##############
        while [[ $temp1 -lt ${p} ]]
        do
          let ts2=${ts2}+${n_states_neut_sym[$temp1]} 
          let temp1=${temp1}+1
        done
           ##############
        while [[ $n -le ${n_states_neut_sym[$p]} ]]
        do
          let ts2=${ts2}+1
          echo "seeking for dipole moment of neutral  ${n}.${pp}-${m}.${op} = ${ts2}-${ts1}"
          a=$(grep -m1 "<${n}.${pp}|DM${LL}|${m}.${op}>" ${fol}/${molpro_input}${R}.out | awk '{print $4}')
          if [[ -z $a ]]; then
            echo "0.0" >> ${output_loc}/LiH_DM${LL}_${ts2}_${ts1}.input
          else
            echo "$a" >> ${output_loc}/LiH_DM${LL}_${ts2}_${ts1}.input
          fi

          let n=$n+1
        done
        let m=$m+1
      done

      fi
      let p=$p+1
    done
    let o=$o+1
  done
done

rm ${fol}/*.cube
mv ${fol}/* ${output_loc}
rm -r ${fol}

let i=$i+1
done

exit

