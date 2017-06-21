#!/usr/bin/env gnuplot



set encoding cp1250 # make the minus sign longer
set terminal push
set terminal postscript eps enhanced font "Times, 22"
set output "pt_hist.eps"
set multiplot


reset

ht1 = 0.29
ht2 = 0.29
ht3 = 1 - ht1 - ht2

set label 1 at screen 0.01, 0.98         "(a)"
set label 2 at screen 0.01, ht2+ht3-0.02 "(b)"
set label 3 at screen 0.01, ht3-0.03     "(c)"


set origin 0, ht2 + ht3
set size 1, ht1

set lmargin 8
set bmargin 0
set tmargin 0.2

set format x ""
set xtics 200
set mxtics 2
set xrange [-1850:-800]
#set xlabel "Shifted simulation time, {/Times-Italic t} + {/Times-Italic t}_{/*0.8 0}"
#
#set format y "10^{/*0.8 %T}"
set ytics 0.01 offset 0.5, 0
set mytics 10
set yrange [0:0.014]
set ylabel "Histogram" offset 2.5, 0


#set key left bottom Left reverse font "Times, 22" spacing 1.5

plot [:][:] \
    "../../data/pt/pt2gaus_L32.his" u 1:($5 >= 0 ? $2 : 1/0) w l lt 1 t "", \
    ""                              u 1:($5 <  0 ? $2 : 1/0) w l lt 2 t "", \
    -1 notitle




set origin 0, ht3
set size 1, ht2

set ytics 0.01
set yrange [1.405:1.445]

set ylabel "{/Times-Italic c}_1/{/Symbol-Oblique s}_{/Times-Italic E}"


plot [:][:] \
    "../../data/pt/pt2gaus_L32.dat" u 2:($4/32) w p pt 4 notitle



set origin 0, 0
set size 1, ht3

set bmargin 2.8

set xtics 200 offset 0, 0.3
set mxtics 2
set format x "%g"
set xlabel "Energy, {/Times-Italic E}" offset 0, 0.3

set ytics 0.1
set yrange [0.55:0.85]

set ylabel "{/Times-Italic c}_2" offset 1.9, 0


plot [:][:] \
    "../../data/pt/pt2gaus_L32.dat" u 2:($5) w p pt 4 notitle



unset multiplot

unset output
set terminal pop
reset
