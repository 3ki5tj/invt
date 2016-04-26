#!/usr/bin/env gnuplot



# Normalized error vs the width of the Gaussian updating scheme



set encoding cp1250 # make the minus sign longer
set terminal push
# dl 2.0 make dashed line longer
set terminal postscript eps enhanced dl 2.0 size 3.5, 4 font "Times, 26"
set output "erriascan.eps"
set multiplot


hbot = 0.45
htop = 1 - hbot

# bottom panel
set size 1, hbot
set origin 0, 0

set border 1+2+8

set tmargin 0.3
set lmargin 7

set logscale x
set xtics offset 0, 0.3 nomirror
set mxtics 10
set format x "10^{/*0.8 %T}"
#set xlabel "{/Symbol-Oblique s} (Gaussian) or {/Times-Italic n}/{/Times-Italic K}/{/Symbol \326}(2{/Symbol-Italic p}) (Bandpass)" offset 0, 0.5
set xlabel "Initial updating magnitude, {/Symbol-Oblique a}(0)" offset 0, 0.5
set xrange [1e-6:2e-2]

#set logscale y
#set format y "10^{/*0.8 %T}"
set ytics 5e-7 offset 0.3, 0
set mytics 10
set format y "%.0t{/Symbol \264}10^{/*0.8 %T}"
#set format y "%.0t{/Symbol \327}10^{%L}"
#set mytics 2
set ylabel "Error, {/Times-Italic E}" offset 1.5, 3.0

a0 = 0.0001
fac = 2*100/sqrt(2*pi)

plot [:][2e-7:1.4e-6] \
    "../../data/iascan/iarun_sig10_g.dat" u 1:($2)    w p pt 13 ps 1.4      notitle, \
    "../../data/iascan/iaprd_sig10_g.dat" u 1:($2**2) w l lt 1 lw 6         notitle, \
    -1 notitle




# top panel
set origin 0, hbot
set size 1, hbot

set tmargin 5
set bmargin 0

# `width` to reduce the text length
set key Left reverse width -7.5 spacing 1.0 font "Times, 20"

unset xlabel
set format x ""
unset xtics

set border 2+4+8
set logscale x2
set x2tics nomirror
set mx2tics 10
set format x2 ""
set x2range [1e-6:2e-2]

set ytics 1e-4
unset ylabel

plot [:][5e-5:3e-4] \
    "../../data/iascan/iarun_sig10_l.dat" u 1:($2)    w p pt 14 ps 1.4   axes x2y1   notitle, \
    "../../data/iascan/iaprd_sig10_l.dat" u 1:($2**2) w l lt 2 lw 6      axes x2y1   notitle, \
    -1 notitle

unset multiplot
unset output
set terminal pop
reset
