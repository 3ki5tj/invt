#!/usr/bin/env gnuplot



# Normalized error vs the proportionality constant c
# as in the schedule of updating magnitude
# alpha(t) = c / (t + t0)
# for the Gaussian updating scheme



set encoding cp1250 # make the minus sign longer
set terminal push
# dl 3.0 make dashed line longer
set terminal postscript eps enhanced dl 3.0 size 5, 3.5 font "Times, 24"
set output "errsig5.eps"


set xtics 5 offset 0, 0.3
set mxtics 5
set xlabel "{/Times 1 / }{/Symbol-Oblique l}" offset 0, 0.5

set logscale y
set format y "10^{/*0.8 %T}"
set mytics 10
set ylabel "Normalized error, {/Times (}{/Times-Italic t}{/Times + }{/Times-Italic t}_{/Times 0} {/Times )} {/Times-Italic E}"

# `width` to reduce the text length
set key at 30, 4e4 Left reverse width -7 spacing 1.4

a0 = 0.0001

plot [0:30][40:1e5] \
    "../../data/invt/sig5/sig5_g_err.dat"     u 1:($2*(1e8+$1/a0))     every 2 w p pt 5 ps 2.0  lc rgb "#000000"    t "Global, {/Times-Italic t} = 10^8, simulation", \
    "../../data/invt/sig5/sig5_g_prd.dat"     u 1:($2**2*(1e8+$1/a0))          w l lt 1 lw 2    lc rgb "#000000"    t "Global, {/Times-Italic t} = 10^8, theory", \
    "../../data/invt/sig5/sig5t1e7_g_err.dat" u 1:($2*(1e7+$1/a0))     every 2 w p pt 5 ps 2.0  lc rgb "#808080"    t "Global, {/Times-Italic t} = 10^7, simulation", \
    "../../data/invt/sig5/sig5t1e7_g_prd.dat" u 1:($2**2*(1e7+$1/a0))          w l lt 1 lw 2    lc rgb "#808080"    t "Global, {/Times-Italic t} = 10^7, theory", \
    "../../data/invt/sig5/sig5_l_err.dat"     u 1:($2*(1e8+$1/a0))     every 2 w p pt 7 ps 2.0  lc rgb "#000000"    t "Local, {/Times-Italic t} = 10^8, simulation", \
    "../../data/invt/sig5/sig5_l_prd.dat"     u 1:($2**2*(1e8+$1/a0))          w l lt 2 lw 2    lc rgb "#000000"    t "Local, {/Times-Italic t} = 10^8, theory", \
    "../../data/invt/sig5/sig5t1e7_l_err.dat" u 1:($2*(1e7+$1/a0))     every 2 w p pt 7 ps 2.0  lc rgb "#808080"    t "Local, {/Times-Italic t} = 10^7, simulation", \
    "../../data/invt/sig5/sig5t1e7_l_prd.dat" u 1:($2**2*(1e7+$1/a0))          w l lt 2 lw 2    lc rgb "#808080"    t "Local, {/Times-Italic t} = 10^7, theory", \
    -1 notitle


unset output
set terminal pop
reset
