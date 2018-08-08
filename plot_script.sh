#!/bin/bash

if [[ -z $1 ]]
then
   file=Output.log
else
   file=$1
fi
grep "Electric field (Z)" $file | sed "s/Electric field (Z)//" > gp_efield.txt
for i in 2 3 4 5 6 7 8 9 10
do
grep "Population on neutral state $i = " $file | sed "s/Population on neutral state $i = //" > gp_pop_$i.txt
done
for i in 1
do
grep "Population on cation state $i = " $file | sed "s/Population on cation state $i = //" > gp_pop_cat_$i.txt
done
grep "Time " $file | sed "s/Time //" > time.txt

gnf=gnufile.gp
lw1=4
lw2=4
lw3=4
lw4=4
lw5=4
lw6=4
lw7=4
lw8=4
lw9=4
lw10=4
lw11=4
lw12=4
color1=black
color2=red
color3=blue
color4=green
color5=purple
color6=gold
color7=magenta
color8=dark-green
color9=pink
color10=violet
color11=dark-blue
color12=black
lt1=0
lt2=1
lt3=1
lt4=1
lt5=1
lt6=1
lt7=1
lt8=1
lt9=1
lt10=1
lt11=1
lt12=2

dt1=2
dt2=1
dt3=1
dt4=1
dt5=1
dt6=1
dt7=1
dt8=1
dt9=1
dt10=1
dt11=1
dt12=1

#set terminal x11 enhanced font "Times, 20 pts"
#set title "populations and e-field"
cat > $gnf << MAFG
set terminal postscript enhanced color font "Times-Roman, 20 pts"
set output "populations_0_10fs.eps"
set xlabel "Time (fs)"
set ylabel "|c_i(t)|^2"
set ylabel offset 1,0
set xrange[0:20]

#set for [i=1:7] linetype i dt i
set style line 1 lt $lt1 lc rgb "${color1}" lw $lw1 dt $dt1 
set style line 2 lt $lt2 lc rgb "${color2}" lw $lw2 dt $dt2
set style line 3 lt $lt3 lc rgb "${color3}" lw $lw3 dt $dt3
set style line 4 lt $lt4 lc rgb "${color4}" lw $lw4 dt $dt4
set style line 5 lt $lt5 lc rgb "${color5}" lw $lw5 dt $dt5
set style line 6 lt $lt6 lc rgb "${color6}" lw $lw6 dt $dt6
set style line 7 lt $lt7 lc rgb "${color7}" lw $lw7 dt $dt7
set style line 8 lt $lt8 lc rgb "${color8}" lw $lw8 dt $dt8
set style line 9 lt $lt9 lc rgb "${color9}" lw $lw9 dt $dt9
set style line 10 lt $lt10 lc rgb "${color10}" lw $lw10 dt $dt10
set style line 11 lt $lt11 lc rgb "${color11}" lw $lw11 dt $dt11
set style line 12 lt $lt12 lc rgb "${color12}" lw $lw12 dt $dt12

set key right top
plot "<paste time.txt gp_efield.txt" u 1:($2)*10 t 'Pulse' w l ls 1 \\
MAFG

for i in 2 3 4 5 6 7 8 9 10 
do
   let j=$i+1
   let m=$i-1
cat >> $gnf << MAFG
,"<paste time.txt gp_pop_$i.txt" u 1:2 t '$m {/Symbol S}' w l ls $i \\
MAFG
done

cat >> $gnf << MAFG
,"<paste time.txt gp_pop_cat_1.txt" u 1:2 t 'Cation' w l ls 12 
MAFG

gnuplot -p "$gnf"

