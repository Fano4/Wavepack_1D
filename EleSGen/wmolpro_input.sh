
function prefiles_input () {
   echo "

***,lithium_hydride
 
gprint,orbitals,civector
geometry={
 2
 LiH
  H          0.0000000000        0.0000000000      ${RH}
  Li         0.0000000000        0.0000000000      ${RLi}
 }
 basis=${basis_set}


file,2,${molpro_wfu_file}${R}.wfu,new

{casscf
maxiter, 40
closed, 0 , 0 , 0 , 0
occ, ${n_occ_sym[0]} , ${n_occ_sym[1]} , ${n_occ_sym[2]} , ${n_occ_sym[3]}
!pspace,5
config
wf,4,1,0;state, ${n_states_neut_sym[0]};weight, ${weight_neut_states_sym[0]};
wf,4,2,0;state, ${n_states_neut_sym[1]};weight, ${weight_neut_states_sym[1]};
wf,4,3,0;state, ${n_states_neut_sym[2]};weight, ${weight_neut_states_sym[2]};
wf,4,4,0;state, ${n_states_neut_sym[3]};weight, ${weight_neut_states_sym[3]};
}

file,2,${molpro_wfu_prefile}${R}.wfu,new

{casscf
maxiter,40
closed, 0
occ, ${n_occ_sym[0]} , ${n_occ_sym[1]} , ${n_occ_sym[2]} , ${n_occ_sym[3]}
!pspace,5
config
wf,3,1,1,1;state, ${n_states_cat_sym[0]} ;
wf,3,2,1,1;state, ${n_states_cat_sym[1]} ;
wf,3,3,1,1;state, ${n_states_cat_sym[2]} ;
!wf,3,4,1,1;state, ${n_states_cat} ;
}

---
   " > ${fol}/${molpro_prefile}.com
}

function main_input () {
   echo "

***,lithium_hydride
 
gprint,civector,orbitals
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

file,2,${molpro_wfu_file}${R}.wfu

{casscf
maxiter, 40
closed, 0 , 0 , 0 , 0
occ, ${n_occ_sym[0]} , ${n_occ_sym[1]} , ${n_occ_sym[2]} , ${n_occ_sym[3]}
!pspace,5
config
wf,4,1,0;state, ${n_states_neut_sym[0]};weight, ${weight_neut_states_sym[0]};
wf,4,2,0;state, ${n_states_neut_sym[1]};weight, ${weight_neut_states_sym[1]};
wf,4,3,0;state, ${n_states_neut_sym[2]};weight, ${weight_neut_states_sym[2]};
wf,4,4,0;state, ${n_states_neut_sym[3]};weight, ${weight_neut_states_sym[3]};
}

!{cube,LiH_neut.cub,-1,125,125,125;
! orbital,occ
! BRAGG,10
! }

file,2,${molpro_wfu_prefile}${R}.wfu

{casscf
maxiter,40
closed, 0
occ, ${n_occ_sym[0]} , ${n_occ_sym[1]} , ${n_occ_sym[2]} , ${n_occ_sym[3]}
!pspace,5
config
wf,3,1,1,1;state, ${n_states_cat_sym[0]} ;
wf,3,2,1,1;state, ${n_states_cat_sym[1]} ;
wf,3,3,1,1;state, ${n_states_cat_sym[2]} ;
!wf,3,4,1,1;state, ${n_states_cat} ;
}

!{casscf
!maxiter,40
!closed, 0
!occ, ${n_occ_sym[0]} , ${n_occ_sym[1]} , ${n_occ_sym[2]} , ${n_occ_sym[3]}
!!pspace,5
!config
!wf,2,1,0,2;state, ${n_states_dicat_sym1} ;
!}
---
   " > ${fol}/${molpro_input}${R}.com
}

function NAC_input () {
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

file,2,${molpro_wfu_file}${R}.wfu

   " > ${fol}/${molpro_nac_input}${R}.com
let sym=0
while [[ $sym -lt ${n_sym}  ]]
do
   let symp=$sym+1
   let m=1
   echo "$symp / ${n_sym} => ${n_states_neut_sym[${sym}]}"
   while [[ $m -le ${n_states_neut_sym[${sym}]} ]]
   do
      let n=1
      while [[ $n -lt $m ]]
      do
         echo "
{casscf
maxiter, 40
closed, 0
occ, ${n_occ_sym[0]} , ${n_occ_sym[1]} , ${n_occ_sym[2]} , ${n_occ_sym[3]}
!pspace,5
config
wf,4,1,0;state, ${n_states_neut_sym[0]};weight, ${weight_neut_states_sym[0]};
wf,4,2,0;state, ${n_states_neut_sym[1]};weight, ${weight_neut_states_sym[1]};
wf,4,3,0;state, ${n_states_neut_sym[2]};weight, ${weight_neut_states_sym[2]};
wf,4,4,0;state, ${n_states_neut_sym[3]};weight, ${weight_neut_states_sym[3]};
CPMCSCF,NACM,$m.$symp,$n.$symp,accu=1d-05}
force
         " >> ${fol}/${molpro_nac_input}${R}.com
         echo "CPMCSCF,NACM,$m.$symp,$n.$symp"

         let n=$n+1
      done
      let m=$m+1
   done
   let sym=$sym+1
done

}
