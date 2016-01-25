#!/usr/bin/env gnuplot



# Normalized error vs the proportionality constant c
# as in the schedule of updating magnitude
# alpha(t) = c / (t + t0)
# for the nearest-neighbor updating scheme



set encoding cp1250 # make the minus sign longer
set terminal push
# dl 3.0 make dashed line longer
set terminal postscript eps enhanced dl 3.0 size 5, 3.5 font "Times, 24"
set output "errnnbr2.eps"


set xtics 1 offset 0, 0.3
set mxtics 2
set xlabel "{/Times 1 / }{/Symbol-Oblique l}" offset 0, 0.5

set logscale y
set format y "10^{%T}"
set mytics 10
set ylabel "Normalized error, {/Times (}{/Times-Italic t}{/Times + }{/Times-Italic t}_{/Times 0}{/Times )} {/Times-Italic E}"

# `width` to reduce the text length
set key at 4, 6e3 Left reverse width -7

plot [0:4][80:5e4] \
    "../../data/nb0.2/nb0.2_g_err.dat"     u 1:($2*1e8):($4*1e8)  w e pt 5 ps 1.2  lc rgb "#000000"    t "Global, {/Times-Italic t} = 10^8, simulation", \
    "../../data/nb0.2/nb0.2_g_prd.dat"     u 1:($2**2*1e8)        w l lt 1 lw 2    lc rgb "#000000"    t "Global, {/Times-Italic t} = 10^8, theory", \
    "../../data/nb0.2/nb0.2t1e7_g_err.dat" u 1:($2*1e7):($4*1e7)  w e pt 5 ps 1.2  lc rgb "#a0a0a0"    t "Global, {/Times-Italic t} = 10^7, simulation", \
    "../../data/nb0.2/nb0.2t1e7_g_prd.dat" u 1:($2**2*1e7)        w l lt 1 lw 2    lc rgb "#a0a0a0"    t "Global, {/Times-Italic t} = 10^7, theory", \
    "../../data/nb0.2/nb0.2_l_err.dat"     u 1:($2*1e8):($4*1e8)  w e pt 7 ps 1.5  lc rgb "#000000"    t "Local, {/Times-Italic t} = 10^8, simulation", \
    "../../data/nb0.2/nb0.2_l_prd.dat"     u 1:($2**2*1e8)        w l lt 2 lw 2    lc rgb "#000000"    t "Local, {/Times-Italic t} = 10^8, theory", \
    "../../data/nb0.2/nb0.2t1e7_l_err.dat" u 1:($2*1e7):($4*1e7)  w e pt 7 ps 1.5  lc rgb "#a0a0a0"    t "Local, {/Times-Italic t} = 10^7, simulation", \
    "../../data/nb0.2/nb0.2t1e7_l_prd.dat" u 1:($2**2*1e7)        w l lt 2 lw 2    lc rgb "#a0a0a0"    t "Local, {/Times-Italic t} = 10^7, theory", \
    -1 notitle


unset output
set terminal pop
reset
