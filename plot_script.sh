#!/bin/bash

if [[ -z $1 ]]
then
   file=Output.log
else
   file=$1
fi
grep "Electric field (Z)" $file | sed "s/Electric field (Z)//" > gp_efield.txt
for i in 2 3 4 5 6 7
do
grep "Population on neutral state $i = " $file | sed "s/Population on neutral state $i = //" > gp_pop_$i.txt
done
for i in 1
do
grep "Population on cation state $i = " $file | sed "s/Population on cation state $i = //" > gp_pop_cat_$i.txt
done
grep "Time " $file | sed "s/Time //" > time.txt

gnf=gnufile.gp
lw1=2
lw2=1
lw3=1
lw4=1
lw5=1
lw6=1
lw7=1
lw8=1
lw9=1
color1=black
color2=red
color3=blue
color4=green
color5=purple
color6=gold
color7=magenta
color8=dark-green
color9=black
lt1=0
lt2=1
lt3=1
lt4=1
lt5=1
lt6=1
lt7=1
lt8=1
lt9=2


cat > $gnf << MAFG
set title "populations and e-field"
#set for [i=1:7] linetype i dt i
set style line 1 lt $lt1 lc rgb "${color1}" lw $lw1
set style line 2 lt $lt2 lc rgb "${color2}" lw $lw2
set style line 3 lt $lt3 lc rgb "${color3}" lw $lw3
set style line 4 lt $lt4 lc rgb "${color4}" lw $lw4
set style line 5 lt $lt5 lc rgb "${color5}" lw $lw5
set style line 6 lt $lt6 lc rgb "${color6}" lw $lw6
set style line 7 lt $lt7 lc rgb "${color7}" lw $lw7
set style line 8 lt $lt8 lc rgb "${color8}" lw $lw8
set style line 9 lt $lt9 lc rgb "${color9}" lw $lw9

set key left top
plot "<paste time.txt gp_efield.txt" u 1:2 t 'Pulse' w l ls 1 \\
MAFG

for i in 2 3 4 5 6 7
do
   let j=$i+1
cat >> $gnf << MAFG
,"<paste time.txt gp_pop_$i.txt" u 1:2 t 'State $i' w l ls $j \\
MAFG
done

cat >> $gnf << MAFG
,"<paste time.txt gp_pop_cat_1.txt" u 1:2 t 'Cation' w l ls 9 
MAFG

gnuplot -p "$gnf"

