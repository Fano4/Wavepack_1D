#!/bin/bash

#WRITE THE COORDINATE IN THE OUTPUT FILE




function print_coordinate () {
echo "${R}">>${output_loc}/coordinates.input
}




#GATHER ENERGIES IN MOLPRO OUTPUT




function print_pes_neutral () {
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
   grep -m1 "!MCSCF STATE ${m}.${kp} Energy" ${fol}/${molpro_input}${R}.out | awk '{print $5}' >> ${output_loc}/${pes_neut}${ts}.input
   let m=$m+1
done
let k=$k+1
done
}




function print_pes_cation () {
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
   grep -m2 "!MCSCF STATE ${m}.${kp} Energy" ${fol}/${molpro_input}${R}.out | tail -n1 | awk '{print $5}' >> ${output_loc}/${pes_cat}${ts}.input
   let m=$m+1
done
let k=$k+1
done
}




#GATHER DIPOLES IN MOLPRO OUTPUT



function print_dipole_neutral () {

#loop between the different components of the vector
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
#             if [[ $LL -eq X ]]; then
#               a=$(grep "SA-MC NACME FOR STATES ${n}.${pp} - ${m}.${op}"  ${fol}/${molpro_nac_input}${R}.out)
#               if [[ -z $a ]]; then
#                  let n=$n+1
#                  continue
#               else
#                  var1=$(grep -a4 "SA-MC NACME FOR STATES ${n}.${pp} - ${m}.${op}"  ${fol}/${molpro_nac_input}${R}.out |tail -n1|awk '{print $4}')
#                  var2=$(grep -a5 "SA-MC NACME FOR STATES ${n}.${pp} - ${m}.${op}"  ${fol}/${molpro_nac_input}${R}.out |tail -n1|awk '{print $4}')
#                  res=$(awk "BEGIN {print ( $mLi / ( $mH + $mH ) ) * $var2 - ( $mLi / ( $mLi + $mH ) ) * $var1 }")
#                 echo "$res" >> ${output_loc}/LiH_NAC_${ts2}_${ts1}.input
#               fi
#             fi
             let n=$n+1
          done
          let m=$m+1
        done
     else #if initial sym is not the same as final sym, we are in an off-diagonal block of the matrix. we have to consider all elements
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
}
