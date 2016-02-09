#!/usr/bin/env gnuplot



# Normalized error vs the width of the Gaussian updating scheme



set encoding cp1250 # make the minus sign longer
set terminal push
# dl 3.0 make dashed line longer
set terminal postscript eps enhanced dl 3.0 size 3.5, 5 font "Times, 26"
set output "errsigscan.eps"
set multiplot


hbot = 0.54
htop = 1 - hbot

# bottom panel
set size 1, hbot
set origin 0, 0

set xtics 2 offset 0, 0.3
set mxtics 2
set xlabel "Width of the Gaussian, {/Symbol-Oblique s}" offset 0, 0.5

set logscale y
set format y "10^{%T}"
#set ytics 1e3
#set format y "%.0t{/Symbol \264}10^{%T}"
#set mytics 2
set ylabel "Normalized error, ({/Times-Italic t} + {/Times-Italic t}_{0}) {/Times-Italic E}"

a0 = 0.0001

set title "Local sampling scheme" offset 0, -0.5

plot [:][6000:1e5] \
    "../../data/scan/sigscan_t1e8_l.dat"   u 1:($4**2) w l lt 1 lw 4 lc rgb "#000000" notitle, \
    "../../data/scan/sigscan_t1e8_l.dat"   u 1:($6**2) w l lt 2 lw 4 lc rgb "#000000" notitle, \
    "../../data/scan/sigscan_t1e10_l.dat"  u 1:($4**2) w l lt 4 lw 6 lc rgb "#808080" notitle, \
    "../../data/scan/sigscan_t1e10_l.dat"  u 1:($6**2) w l lt 5 lw 4 lc rgb "#808080" notitle, \
    -1 notitle



# top panel
set origin 0, hbot
set size 1, hbot

set tmargin 4
set bmargin 0

# `width` to reduce the text length
set key Left reverse width -7 spacing 1.0 font "Times, 20"

unset xlabel
set format x ""

set title "Global sampling scheme"

plot [:][30:2e4] \
    "../../data/scan/sigscan_t1e8_g.dat"   u 1:($4**2) w l lt 1 lw 4 lc rgb "#000000"    t "{/Times-Italic t} = 10^{8}, Inv. {/Times t}", \
    "../../data/scan/sigscan_t1e8_g.dat"   u 1:($6**2) w l lt 2 lw 4 lc rgb "#000000"    t "{/Times-Italic t} = 10^{8}, Optimal", \
    "../../data/scan/sigscan_t1e10_g.dat"  u 1:($4**2) w l lt 4 lw 6 lc rgb "#808080"    t "{/Times-Italic t} = 10^{10}, Inv. {/Times t}", \
    "../../data/scan/sigscan_t1e10_g.dat"  u 1:($6**2) w l lt 5 lw 4 lc rgb "#808080"    t "{/Times-Italic t} = 10^{10}, Optimal", \
    -1 notitle


unset multiplot
unset output
set terminal pop
reset
