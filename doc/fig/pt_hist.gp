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


set style line 1 lt 1
set style line 2 lt 2 lw 3
set style line 3 lt 3 lw 1

set origin 0, ht2 + ht3
set size 1, ht1

set lmargin 8
set bmargin 0
set tmargin 0.2
set rmargin 2.5

set format x ""
set xtics 1000
set mxtics 5
set xrange [-7600:-3000]
#set xlabel "Shifted simulation time, {/Times-Italic t} + {/Times-Italic t}_{/*0.8 0}"
#
#set format y "10^{/*0.8 %T}"
set ytics 0.002 offset 0.5, 0
set mytics 2
set yrange [0:0.007]
set ylabel "Histogram" offset 2.5, 0


#set key left bottom Left reverse font "Times, 22" spacing 1.5

plot [:][:] \
    "../../data/pt/pt2gaus_L64b.his" u 1:($5 >= 0 ? $2 : 1/0) w l ls 1 t "", \
    "../../data/pt/lng_L64b.dat"     u 1:($4)                 w l ls 2 t "", \
    -1 notitle




set origin 0, ht3
set size 1, ht2

set ytics 0.01
set yrange [1.405:1.445]

set ylabel "{/Times-Italic c}_1/{/Symbol-Oblique s}_{/Times-Italic E}"


plot [:][:] \
    "../../data/pt/pt2gaus_L64b.dat" u 2:($4/64) w p pt 4 lc "black" notitle


#    "../../data/pt/lng_L64b.dat"     u 1:($1 < -3200 ? $3 : 1/0)                 w l ls 2 lw 1 t "", \


set origin 0, 0
set size 1, ht3

set bmargin 2.8

set xtics 1000 offset 0, 0.3
set mxtics 5
set format x "%g"
set xlabel "Energy, {/Times-Italic E}" offset 0, 0.3

set ytics 0.5
set mytics 5
set yrange [0.3:1.3]

set ylabel "{/Symbol \326}~2{.7-} {/Times-Italic c}_2" offset 1.9, 0


plot [:][:] \
    1 w lp lt -1 pi -2 lw 1.0 lc rgb "#888888" notitle, \
    "../../data/pt/pt2gaus_L64b.dat" u 2:($5*sqrt(2)) w p pt 4 lc "black" notitle



unset multiplot

unset output
set terminal pop
reset
