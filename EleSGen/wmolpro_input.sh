
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


file,2,${molpro_wfu_file}${R}.wfu
restart,2;
geometry={
 2
 LiH
  H          0.0000000000        0.0000000000      ${RH}
  Li         0.0000000000        0.0000000000      ${RLi}
 }

{casscf
maxiter, 40
closed, 0 , 0 , 0 , 0
occ, ${n_occ_sym[0]} , ${n_occ_sym[1]} , ${n_occ_sym[2]} , ${n_occ_sym[3]}
!pspace,5
config
wf,4,1,0;state, ${n_states_neut_sym[0]};!weight, ${weight_neut_states_sym[0]};
wf,4,2,0;state, ${n_states_neut_sym[1]};!weight, ${weight_neut_states_sym[1]};
wf,4,3,0;state, ${n_states_neut_sym[2]};!weight, ${weight_neut_states_sym[2]};
wf,4,4,0;state, ${n_states_neut_sym[3]};!weight, ${weight_neut_states_sym[3]};
ORBITAL,2140.2
}

file,2,${molpro_wfu_prefile}${R}.wfu
restart,2;
geometry={
 2
 LiH
  H          0.0000000000        0.0000000000      ${RH}
  Li         0.0000000000        0.0000000000      ${RLi}
 }

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
 
gprint,basis,civector,orbitals
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
restart,2;

{casscf
maxiter, 40
closed, 0 , 0 , 0 , 0
occ, ${n_occ_sym[0]} , ${n_occ_sym[1]} , ${n_occ_sym[2]} , ${n_occ_sym[3]}
!pspace,5
config
wf,4,1,0;state, ${n_states_neut_sym[0]};!weight, ${weight_neut_states_sym[0]};
wf,4,2,0;state, ${n_states_neut_sym[1]};!weight, ${weight_neut_states_sym[1]};
wf,4,3,0;state, ${n_states_neut_sym[2]};!weight, ${weight_neut_states_sym[2]};
wf,4,4,0;state, ${n_states_neut_sym[3]};!weight, ${weight_neut_states_sym[3]};
ORBITAL,2140.2
}

 {matrop;
    load,mo,orb,2140.2;
    load,dmx,oper,dmx
    tran,dmx_mo,dmx,mo
    print dmx_mo
 }
 {matrop;
    load,mo,orb,2140.2;
    load,dmy,oper,dmy
    tran,dmy_mo,dmy,mo
    print dmy_mo
 }
 {matrop;
    load,mo,orb,2140.2;
    load,dmz,oper,dmz
    tran,dmz_mo,dmz,mo
    print dmz_mo
 }

file,2,${molpro_wfu_prefile}${R}.wfu
restart,2;

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

***, Lithium_hydride_NAC

R=${R}
dR=${dR}
mH=${mH}
mLi=${mLi}

RH=-(R)*(mLi/(mH+mLi))
RLi=(R)*(mH/(mH+mLi))

geometry={
 2
 LiH
  H          0.0000000000        0.0000000000      RH
  Li         0.0000000000        0.0000000000      RLi
 }

 basis=${basis_set}

   " > ${fol}/${molpro_nac_input}${R}.com

let sym=0
while [[ $sym -lt ${n_sym}  ]]
do
   let symp=$sym+1

echo " 
file,2,${molpro_wfu_file}${R}.wfu

{casscf
maxiter, 40
closed, 0 , 0 , 0 , 0
occ, ${n_occ_sym[0]} , ${n_occ_sym[1]} , ${n_occ_sym[2]} , ${n_occ_sym[3]}
!pspace,5
config
wf,4,1,0;state, ${n_states_neut_sym[0]};!weight, ${weight_neut_states_sym[0]};
wf,4,2,0;state, ${n_states_neut_sym[1]};!weight, ${weight_neut_states_sym[1]};
wf,4,3,0;state, ${n_states_neut_sym[2]};!weight, ${weight_neut_states_sym[2]};
wf,4,4,0;state, ${n_states_neut_sym[3]};!weight, ${weight_neut_states_sym[3]};
orbital,2140.2
}

{ci
wf,4,${symp},0;state, ${n_states_neut_sym[${sym}]};
noexc
save,6000.2
dm,8${sym}00.2
}

RH=-(R+dR)*(mLi/(mH+mLi))
RLi=(R+dR)*(mH/(mH+mLi))
geometry={
 2
 LiH
  H          0.0000000000        0.0000000000      RH
  Li         0.0000000000        0.0000000000      RLi
 }
{casscf
maxiter, 40
closed, 0 , 0 , 0 , 0
occ, ${n_occ_sym[0]} , ${n_occ_sym[1]} , ${n_occ_sym[2]} , ${n_occ_sym[3]}
!pspace,5
config
wf,4,1,0;state, ${n_states_neut_sym[0]};!weight, ${weight_neut_states_sym[0]};
wf,4,2,0;state, ${n_states_neut_sym[1]};!weight, ${weight_neut_states_sym[1]};
wf,4,3,0;state, ${n_states_neut_sym[2]};!weight, ${weight_neut_states_sym[2]};
wf,4,4,0;state, ${n_states_neut_sym[3]};!weight, ${weight_neut_states_sym[3]};
orbital,2141.2
diab,2140.2
}
{ci
wf,4,${symp},0;state, ${n_states_neut_sym[${sym}]};
noexc
save,6001.2
}
{ci;trans,6000.2,6001.2;
dm,8${sym}01.2}

RH=-(R-dR)*(mLi/(mH+mLi))
RLi=(R-dR)*(mH/(mH+mLi))
geometry={
 2
 LiH
  H          0.0000000000        0.0000000000      RH
  Li         0.0000000000        0.0000000000      RLi
 }
{casscf
maxiter, 40
closed, 0 , 0 , 0 , 0
occ, ${n_occ_sym[0]} , ${n_occ_sym[1]} , ${n_occ_sym[2]} , ${n_occ_sym[3]}
!pspace,5
config
wf,4,1,0;state, ${n_states_neut_sym[0]};!weight, ${weight_neut_states_sym[0]};
wf,4,2,0;state, ${n_states_neut_sym[1]};!weight, ${weight_neut_states_sym[1]};
wf,4,3,0;state, ${n_states_neut_sym[2]};!weight, ${weight_neut_states_sym[2]};
wf,4,4,0;state, ${n_states_neut_sym[3]};!weight, ${weight_neut_states_sym[3]};
orbital,2142.2
diab,2140.2
}

{ci
wf,4,${symp},0;state, ${n_states_neut_sym[${sym}]};
noexc
save,6002.2
}

{ci;trans,6000.2,6002.2;
dm,8${sym}02.2}
" >> ${fol}/${molpro_nac_input}${R}.com
   let m=1
   while [[ $m -le ${n_states_neut_sym[${sym}]} ]]
   do
      let n=1
      while [[ $n -lt $m ]]
      do
echo "
{ddr,2*dR
state,${m}.${symp},${n}
orbital,2140.2,2141.2,2142.2;
density,8${sym}00.2,8${sym}01.2,8${sym}02.2}
   " >> ${fol}/${molpro_nac_input}${R}.com
         let n=$n+1
      done
      let m=$m+1
   done

   let sym=$sym+1
done

}

function out_and_cube_gen () {

echo "
***, LiH_NAC_computation

file,2,${molpro_wfu_file}${1}.wfu
{cube,${2},-1,125,125,125;
orbital,occ
BRAGG,10}
" >> ${3} 

}
function overlap_geom () {

   R1H=$(awk " BEGIN { print - ( $mLi / ( $mH + $mLi ) ) * $1 } " )
   R1Li=$(awk " BEGIN { print ( $mH / ( $mH + $mLi ) ) * $1 } " )
   R2H=$(awk " BEGIN { print - ( $mLi / ( $mH + $mLi ) ) * $R } " )
   R2Li=$(awk " BEGIN { print ( $mH / ( $mH + $mLi ) ) * $R } " )

   echo "***,molpro_ovlp_computation
   gprint,basis
   GTHRESH,THROVL=-5
   geometry={
   4
   LiH
   H    0.0000    0.0000    ${R1H}
   Li   0.0000    0.0000    ${R1Li}
   H    0.0000    0.0000    ${R2H}
   Li   0.0000    0.0000    ${R2Li}
   }
   basis=${basis_set}
   {matrop;
   load,s
   print s
   }
   ---">$3

}
