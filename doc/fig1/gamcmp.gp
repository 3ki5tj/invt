#!/usr/bin/env gnuplot



# comparison of the Gamma values



set encoding cp1250 # make the minus sign longer
set terminal push
# dl 4.0 make dashed line longer
set terminal postscript eps enhanced dl 4.0 size 3.5, 2.5 font "Times, 24"
set output "gamcmp.eps"

set xtics 20 offset 0, 0.3
set mxtics 4
set xlabel "Eigenmode" offset 0, 0.5

#set logscale x
#set format x "10^{/*0.8 %T}"
#set xtics offset 0, 0.2
#set xlabel "{/Symbol G}^{&{i}one{/*0.7 -}step}"
#set xlabel "{/Symbol G}^{&{i}Local MC}"

set logscale y
set format y "10^{/*0.8 %T}"
set ylabel "{/Symbol G} / {/Symbol G}_{/*0.8 max}"

#set key Left reverse width -2 spacing 1.8

fngamma = "../../data/gamcmp/gamma.dat"

#f(x) = a * x + b
#fit f(x) fngamma u 2:4 via a, b

max2 = `cut -f1 ../../data/gamcmp/gamma.dat | sed 's/  */\t/g' | cut -f3 | head -n 1`
max4 = `cut -f2 ../../data/gamcmp/gamma.dat | sed 's/  */\t/g' | cut -f3 | head -n 1`
max6 = `cut -f3 ../../data/gamcmp/gamma.dat | sed 's/  */\t/g' | cut -f3 | head -n 1`

set key spacing 1.5

plot [:][1e-6:] \
  fngamma u ($1):($2/max2) w p pt 1 ps 0.8 t "One{/*0.7 -}step", \
  fngamma u ($1):($6/max6) w p pt 6 ps 0.8 t "Local MC", \
  fngamma u ($1):($4/max4) w p pt 8 ps 0.8 t "MD"


unset output
set terminal pop
reset
