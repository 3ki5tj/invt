#!/usr/bin/env gnuplot



# Normalized error vs the width of the Gaussian updating scheme



set encoding cp1250 # make the minus sign longer
set terminal push
# dl 3.0 make dashed line longer
set terminal postscript eps enhanced dl 3.0 size 3.5, 5 font "Times, 26"
set output "errsigscan.eps"
set multiplot


hbot = 0.56
htop = 1 - hbot

# bottom panel
set size 1, hbot
set origin 0, 0

set xtics 2 offset 0, 0.3
set mxtics 2
set xlabel "Width of the Gaussian, {/Symbol-Oblique s}" offset 0, 0.5
set xrange [0:12]

set logscale y
set format y "10^{%T}"
#set ytics 1e3
#set format y "%.0t{/Symbol \264}10^{%T}"
#set mytics 2
set ylabel "Normalized error, ({/Times-Italic T} + {/Times-Italic t}_{0}) {/Times-Italic E}"

a0 = 0.0001

set title "MD" offset 0, -0.5

plot [:][5e4:1e6] \
    "../../data/scan/sigscan_t1e8_md.dat"  u 1:($6**2) w l lt 1 lw 2 lc rgb "#000000" notitle, \
    "../../data/scan/sigscan_t1e8_md.dat"  u 1:($4**2) w l lt 2 lw 2 lc rgb "#000000" notitle, \
    "../../data/scan/sigscan_t1e10_md.dat" u 1:($6**2) w l lt 1 lw 6 lc rgb "#000000" notitle, \
    "../../data/scan/sigscan_t1e10_md.dat" u 1:($4**2) w l lt 2 lw 6 lc rgb "#000000" notitle, \
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

set title "MC, global"

plot [:][30:3e4] \
    "../../data/scan/sigscan_t1e8_g.dat"   u 1:($6**2) w l lt 1 lw 2 lc rgb "#000000"    t "{/Times-Italic T} = 10^{8}, optimal", \
    "../../data/scan/sigscan_t1e8_g.dat"   u 1:($4**2) w l lt 2 lw 2 lc rgb "#000000"    t "{/Times-Italic T} = 10^{8}, inv. {/Times t}", \
    "../../data/scan/sigscan_t1e10_g.dat"  u 1:($6**2) w l lt 1 lw 6 lc rgb "#000000"    t "{/Times-Italic T} = 10^{10}, optimal", \
    "../../data/scan/sigscan_t1e10_g.dat"  u 1:($4**2) w l lt 2 lw 6 lc rgb "#000000"    t "{/Times-Italic T} = 10^{10}, inv. {/Times t}", \
    -1 notitle


unset multiplot
unset output
set terminal pop
reset
