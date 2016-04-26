#!/usr/bin/env gnuplot



# Normalized error vs the width of the Gaussian updating scheme



set encoding cp1250 # make the minus sign longer
set terminal push
# dl 2.0 make dashed line longer
set terminal postscript eps enhanced dl 2.0 size 3.5, 2.5 font "Times, 26"
set output "erriascan.eps"
set multiplot

set logscale x
set xtics offset 0, 0.3 nomirror
set mxtics 10
set format x "10^{/*0.8 %T}"
set xlabel "Initial updating magnitude, {/Symbol-Oblique a}(0)" offset 0, 0.5
set xrange [1e-6:1e-2]

set ytics 0.5 offset 0.3, 0
set mytics 5
set ylabel "Relative error, {/Times-Italic E} / {/Times-Italic E}{/*0.6 min}" offset 2.5, 0.0

set key at 5e-3, 1.95 Left reverse

a0 = 0.0001
fac = 2*100/sqrt(2*pi)

emin_g = 0.0005192 ** 2 
emin_l = 0.009561 ** 2

plot [:][0.8:2] \
    "../../data/iascan/iarun_sig10_g.dat" u 1:($2    / emin_g)  w p pt 13 ps 1.4      notitle, \
    "../../data/iascan/iaprd_sig10_g.dat" u 1:($2**2 / emin_g)  w l lt 1 lw 4         t "Global", \
    "../../data/iascan/iarun_sig10_l.dat" u 1:($2    / emin_l)  w p pt 12 ps 1.4      notitle, \
    "../../data/iascan/iaprd_sig10_l.dat" u 1:($2**2 / emin_l)  w l lt 2 lw 4         t "Local", \
    -1 notitle


# phantom plot for the legend
set key at 2e-4, 1.95
plot [:][0.8:2] \
    -1 w p pt 13 ps 1.4 t " ", \
    -1 w p pt 12 ps 1.4 t " "



unset multiplot
unset output
set terminal pop
reset
