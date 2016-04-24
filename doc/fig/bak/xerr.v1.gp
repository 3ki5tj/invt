#!/usr/bin/env gnuplot



# Normalized error vs the width of the Gaussian updating scheme



set encoding cp1250 # make the minus sign longer
set terminal push
# dl 2.0 make dashed line longer
set terminal postscript eps enhanced dl 2.0 size 3.5, 5 font "Times, 26"
set output "xerr.eps"
set multiplot


hbot = 0.56
htop = 1 - hbot

# bottom panel
set size 1, hbot
set origin 0, 0

set xtics 5 offset 0, 0.3
set mxtics 5
#set xlabel "{/Symbol-Oblique s} (Gaussian) or {/Times-Italic n}/{/Times-Italic K}/{/Symbol \326}(2{/Symbol-Italic p}) (Bandpass)" offset 0, 0.5
set xlabel "Eigenmode" offset 0, 0.5
set xrange [0:30]

set logscale y
set format y "10^{/*0.7 %T}"
set ytics 10000
#set format y "%.0t{/Symbol \264}10^{%T}"
#set mytics 2
set ylabel "Component of error"

a0 = 0.0001
fac = 2*100/sqrt(2*pi)

set title "Local" offset 0, -0.5

set key Left reverse spacing 1.1 width -1

plot [:][1e-14:] \
    "../../data/xerr/xerr_l.dat"  u 1:($2) w lp lt 1 lw 1 pt 6  t "Init.", \
    "../../data/xerr/xerr_l.dat"  u 1:($3) w lp lt 2 lw 2 pt 7  t "Final", \
    "../../data/xerr/xerr_l.dat"  u 1:($4) w lp lt 3 lw 2 pt 11 t "Res.", \
    "../../data/xerr/xerr_l.dat"  u 1:($5) w lp lt 4 lw 2 pt 9  t "Asym.", \
    -1 notitle



# top panel
set origin 0, hbot
set size 1, hbot

set tmargin 5
set bmargin 0

unset xlabel
set format x ""

set title "Global"



plot [:][1e-18:1e-3] \
    "../../data/xerr/xerr_g.dat"  u 1:($2) w lp lt 1 lw 1 pt 6  t "Init.", \
    "../../data/xerr/xerr_g.dat"  u 1:($3) w lp lt 2 lw 2 pt 7  t "Final", \
    "../../data/xerr/xerr_g.dat"  u 1:($4) w lp lt 3 lw 2 pt 11 t "Res.", \
    "../../data/xerr/xerr_g.dat"  u 1:($5) w lp lt 4 lw 2 pt 9  t "Asym.", \
    -1 notitle



unset multiplot
unset output
set terminal pop
reset
