#!/bin/bash

function nac_findiff_gen () {

   tempm2=tempm2.txt
   tempm1=tempm1.txt
   tempR=tempR.txt
   tempp1=tempp1.txt
   tempp2=tempp2.txt
   temptot=temptot.txt
   nactau=nac_tau.txt
   nacg=nac_g.txt

   bau=0.529 # Bohr to Angstrom

   molprom2=${fol}/molprom2
   molprom1=${fol}/molprom1
   molproR=${fol}/molproR
   molprop1=${fol}/molprop1
   molprop2=${fol}/molprop2

   Rm2d=$(awk " BEGIN { print $Rmin + ( $i - 2 ) * ( $Rmax - $Rmin ) / ($grid_size) } " )
   Rm1d=$(awk " BEGIN { print $Rmin + ( $i - 1 ) * ( $Rmax - $Rmin ) / ($grid_size) } " )
   Rp1d=$(awk " BEGIN { print $Rmin + ( $i + 1 ) * ( $Rmax - $Rmin ) / ($grid_size) } " )
   Rp2d=$(awk " BEGIN { print $Rmin + ( $i + 2 ) * ( $Rmax - $Rmin ) / ($grid_size) } " )
   R=$(awk " BEGIN { print $Rmin + $i * ( $Rmax - $Rmin ) / ($grid_size) } " )

   let n_states_neut_tot=${n_states_neut_sym[0]}+${n_states_neut_sym[1]}+${n_states_neut_sym[2]}+${n_states_neut_sym[3]}

   if [[ $i -ge 2 ]];then
      if [[ $i -le ${grid_size}-3 ]];then
         testv=1
      else
         testv=0
      fi
   else
         testv=0
   fi

   if [[ $testv -eq 1 ]]; then
      echo " i = $i => call NAC "

      #CALL MOLPRO FOR OVERLAP COMPUTATION
      cp ${output_loc}/${molpro_input}${R}.out ${fol}
      cp ${output_loc}/${molpro_input}${Rm2d}.out ${fol}
      cp ${output_loc}/${molpro_input}${Rm1d}.out ${fol}
      cp ${output_loc}/${molpro_input}${Rp1d}.out ${fol}
      cp ${output_loc}/${molpro_input}${Rp2d}.out ${fol}

      overlap_geom ${Rm2d} ${R} ${molprom2}.com
      overlap_geom ${Rm1d} ${R} ${molprom1}.com
      overlap_geom ${Rp1d} ${R} ${molprop1}.com
      overlap_geom ${Rp2d} ${R} ${molprop2}.com

      cd ${fol}
      echo "running molpro..."
      ${molpro_dir} -n4 -m4000M  -d ${fol} -s ${molprom2}.com
      ${molpro_dir} -n4 -m4000M  -d ${fol} -s ${molprom1}.com
      ${molpro_dir} -n4 -m4000M  -d ${fol} -s ${molprop1}.com
      ${molpro_dir} -n4 -m4000M  -d ${fol} -s ${molprop2}.com
      echo "molpro done"

      echo "running overlap in overlap mode..."
      ${overlap_code_dir}/wf_overlap.exe n_sym ${n_sym} n_elec ${n_elec} n_states ${n_states_neut_sym[0]} ${n_states_neut_sym[1]} ${n_states_neut_sym[2]} ${n_states_neut_sym[3]} molpro_output_1 ${molpro_input}${R}.out molpro_output_2 ${molpro_input}${Rm2d}.out  molpro_ovlp ${molprom2}.out > ${fol}/${tempm2}

      ${overlap_code_dir}/wf_overlap.exe n_sym ${n_sym} n_elec ${n_elec} n_states ${n_states_neut_sym[0]} ${n_states_neut_sym[1]} ${n_states_neut_sym[2]} ${n_states_neut_sym[3]} molpro_output_1 ${molpro_input}${R}.out molpro_output_2 ${molpro_input}${Rm1d}.out  molpro_ovlp ${molprom1}.out > ${fol}/${tempm1} 

      ${overlap_code_dir}/wf_overlap.exe n_sym ${n_sym} n_elec ${n_elec} n_states ${n_states_neut_sym[0]} ${n_states_neut_sym[1]} ${n_states_neut_sym[2]} ${n_states_neut_sym[3]} molpro_output_1 ${molpro_input}${R}.out molpro_output_2 ${molpro_input}${Rp1d}.out molpro_ovlp ${molprop1}.out > ${fol}/${tempp1}

      ${overlap_code_dir}/wf_overlap.exe n_sym ${n_sym} n_elec ${n_elec} n_states ${n_states_neut_sym[0]} ${n_states_neut_sym[1]} ${n_states_neut_sym[2]} ${n_states_neut_sym[3]} molpro_output_1 ${molpro_input}${R}.out molpro_output_2 ${molpro_input}${Rp2d}.out molpro_ovlp ${molprop2}.out > ${fol}/${tempp2}
      
      let p=0
      while [[ $p -lt ${n_states_neut_tot} ]]
      do
         let q=0
         while [[ $q -lt ${n_states_neut_tot} ]]
         do
            if [[ $p -eq $q ]]
            then
               echo "1" >> ${fol}/${tempR}
            else 
               echo "0" >> ${fol}/${tempR}
            fi
            let q=$q+1
         done
         let p=$p+1
      done

      cd ${home}

   paste ${fol}/${tempm2} ${fol}/${tempm1} ${fol}/${tempR} ${fol}/${tempp1} ${fol}/${tempp2} ${output_loc}/${phase_output}${R}.txt > ${fol}/${temptot} 
   cat ${fol}/${temptot} | awk " { print 2 * ( (\$1 - 8 * \$2 + 8 * \$4 - \$5) / ( 12 * ( $R - $Rm1d ) ) ) * ${bau}} " > ${fol}/${nactau}
   cat ${fol}/${temptot} | awk " { print ( (-\$1 + 16 * \$2 - 30 * \$3 + 16 * \$4 - \$5) / ( 12 * ( $R - $Rm1d ) ) ) * ${bau} } " > ${fol}/${nacg}
#   cat ${fol}/${temptot} | awk " { print (\$1 - 8 * \$2 + 8 * \$4 - \$5) / ( 12 * ( $R - $Rm1d ) * ${mu}) } "
#   echo "####"
#   cat ${fol}/${temptot} | awk " { print (-\$1 + 16 * \$2 - 30 * \$3 + 16 * \$4 - \$5) / ( 2 * 12 * ( $R - $Rm1d ) * ${mu}) } " 

      let m=0
      while [[ $m -lt ${n_states_neut_tot} ]] 
      do
         let n=$m
         while [[ $n -lt ${n_states_neut_tot} ]]
         do
            let mp=${m}+1
            let np=${n}+1
            let temp1=$(awk "BEGIN {print ${n_states_neut_tot} * ${m} + ${np} }")
            echo "extracting NAC for states ${n} / ${m} => temp1 = ${temp1}"
            head -n ${temp1} ${fol}/${nactau}| tail -n 1
            head -n ${temp1} ${fol}/${nactau}| tail -n 1  >> ${output_loc}/${nac_tau_output}${np}_${mp}.input
            head -n ${temp1} ${fol}/${nacg}| tail -n 1  >> ${output_loc}/${nac_g_output}${np}_${mp}.input

            let n=$n+1
         done
         let m=$m+1
      done

  # rm ${fol}/*.cube ${fol}/*.wfu ${fol}/*.out ${fol}/*.xml ${fol}/.com

else 

      let m=0
      while [[ $m -lt ${n_states_neut_tot} ]] 
      do
         let n=$m
         while [[ $n -lt ${n_states_neut_tot} ]]
         do
            let mp=${m}+1
            let np=${n}+1
            let temp1=$(awk "BEGIN {print ${n_states_neut_tot} * ${m} + ${np} }")
            echo "0.0"  >> ${output_loc}/${nac_tau_output}${np}_${mp}.input
            echo "0.0"  >> ${output_loc}/${nac_g_output}${np}_${mp}.input

            let n=$n+1
         done
         let m=$m+1
      done
   #mv ${fol}/* ${output_loc}
fi

}

function get_phase() {

   Rm1d=$(awk " BEGIN { print $Rmin + ( $i - 1 ) * ( $Rmax - $Rmin ) / ($grid_size) } " )
   R=$(awk " BEGIN { print $Rmin + $i * ( $Rmax - $Rmin ) / ($grid_size) } " )
   molprom1=${fol}/molprom1
   molproR=${fol}/molproR
   let n_states_neut_tot=${n_states_neut_sym[0]}+${n_states_neut_sym[1]}+${n_states_neut_sym[2]}+${n_states_neut_sym[3]}

   if [[ $i -ge 1 ]];then
      echo "getting phase for point $i with R = $R"

      cp ${output_loc}/${molpro_input}${R}.out ${fol}
      cp ${output_loc}/${molpro_input}${Rm1d}.out ${fol}
      overlap_geom ${Rm1d} ${R} ${molprom1}.com

      cd ${fol}
      echo "running molpro..."
      ${molpro_dir} -n4 -m4000M  -d ${fol} -s ${molprom1}.com
      echo "molpro done"

      echo "running overlap in phase mode..."

      ${overlap_code_dir}/wf_overlap.exe n_sym ${n_sym} n_elec ${n_elec} n_states ${n_states_neut_sym[0]} ${n_states_neut_sym[1]} ${n_states_neut_sym[2]} ${n_states_neut_sym[3]} molpro_output_1 ${molpro_input}${R}.out molpro_output_2 ${molpro_input}${Rm1d}.out  molpro_ovlp ${molprom1}.out > ${fol}/temp_phase.txt phase

      if [[ $i -gt 1 ]]; then
         paste ${output_loc}/${phase_output}${Rm1d}.txt ${fol}/temp_phase.txt | awk '{print $1 * $2}' > ${output_loc}/${phase_output}${R}.txt 
      else
         mv ${fol}/temp_phase.txt ${output_loc}/${phase_output}${R}.txt
      fi
      cd ${home}
   else
      let m=0
      while [[ $m -lt ${n_states_neut_tot} ]]
      do
         echo "1" >> ${output_loc}/${phase_output}${R}.txt
         let m=$m+1
      done

   fi


   let m=0
   while [[ $m -lt ${n_states_neut_tot} ]] 
   do
      let mp=${m}+1
      head -n ${mp} ${output_loc}/${phase_output}${R}.txt | tail -n 1  >> ${output_loc}/${phase_output}${mp}.input
      let m=$m+1
   done

#      rm ${fol}/*.cube ${fol}/*.wfu ${fol}/*.out ${fol}/*.xml ${fol}/*.com
#      mv ${fol}/* ${output_loc}
}
