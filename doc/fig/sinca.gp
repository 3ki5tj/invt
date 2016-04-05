#!/usr/bin/env gnuplot



# optimized schedule of the bandpass updating schemes



set encoding cp1250 # make the minus sign longer
set terminal push
# dl 4 make dashed line longer
set terminal postscript eps enhanced dl 4 size 3.5, 2.5 font "Times, 24"
set output "sinca.eps"


set logscale x
set mxtics 10
set format x "10^{%T}"
set xrange [1e4:1e8]
set xlabel "Simulation time, {/Times-Italic t}" offset 0, 0

set logscale y
set format y "10^{%T}"
set mytics 10
#set yrange [1e4:1e9]
set ylabel "{/Symbol-Oblique a} ({/Times-Italic t})" offset 0, 0

set key left bottom Left reverse noinvert spacing 1.5

plot [][2e-9:5e-5] \
    1/(x + 20000) lw 2 t "1 / ({/Times-Italic t} + 2/{/Symbol-Oblique a}_0)", \
    "../../data/scan/sinc_pbc_alpha.dat"    u 1:($2) every 300 w p pt 6 ps 1.6 lw 2 t "Periodic, {/Times-Italic K} = 4", \
    "../../data/scan/sinc_nonpbc_alpha.dat" u 1:($2) every 300 w p pt 2 ps 1.2 lw 2 t "Non{/*0.7 -}periodic, {/Times-Italic K} = 8", \
    -1 notitle



unset multiplot
unset output
set terminal pop
reset
