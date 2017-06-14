#!/usr/bin/env gnuplot



# Normalized error vs the width of the Gaussian updating scheme



set encoding cp1250 # make the minus sign longer
set terminal push
# dl 3.0 make dashed line longer
set terminal postscript eps enhanced dl 3.0 size 8, 4 font "Times, 26"
set output "errsigscan.eps"
set multiplot


wleft = 0.53
wright = 1 - wleft

set size wleft, 1
set origin 0, 0

set xtics 2 offset 0, 0.3
set mxtics 2
set xlabel "Width of the Gaussian, {/Symbol-Oblique s}" offset 0, 0.5 font "Times, 36"

set logscale y
set format y "10^{%T}"
#set ytics 1e3
#set format y "%.0t{/Symbol \264}10^{%T}"
#set mytics 2
set ylabel "Normalized error, ({/Times-Italic t} + {/Times-Italic t}_{0 }) {/Times-Italic E}" font "Times, 36"

# `width` to reduce the text length
set key Left reverse width -6.5 spacing 1.4

set title "Global sampling scheme" font "Times, 40"

a0 = 0.0001

plot [:][30:] \
    "../../data/scan/sigscan_t1e8_g.dat"   u 1:($4**2) w l lt 1 lw 4 lc rgb "#000000"    t "{/Times-Italic t} = 10^{8}, inverse-time", \
    "../../data/scan/sigscan_t1e8_g.dat"   u 1:($6**2) w l lt 2 lw 4 lc rgb "#000000"    t "{/Times-Italic t} = 10^{8}, exact", \
    "../../data/scan/sigscan_t1e10_g.dat"  u 1:($4**2) w l lt 4 lw 6 lc rgb "#808080"    t "{/Times-Italic t} = 10^{10}, inverse-time", \
    "../../data/scan/sigscan_t1e10_g.dat"  u 1:($6**2) w l lt 5 lw 4 lc rgb "#808080"    t "{/Times-Italic t} = 10^{10}, exact", \
    -1 notitle




set origin wleft, 0
set size wright, 1

#set ytics 2e4
#set mytics 5
unset ylabel

# `width` to reduce the text length
set key left Left reverse width -5 spacing 1.4

set title "Local sampling scheme"

set lmargin 2

plot [:][6000:3e5] \
    "../../data/scan/sigscan_t1e8_l.dat"   u 1:($4**2) w l lt 1 lw 4 lc rgb "#000000"    t "{/Times-Italic t} = 10^{8}, inverse-time", \
    "../../data/scan/sigscan_t1e8_l.dat"   u 1:($6**2) w l lt 2 lw 4 lc rgb "#000000"    t "{/Times-Italic t} = 10^{8}, exact", \
    "../../data/scan/sigscan_t1e10_l.dat"  u 1:($4**2) w l lt 4 lw 6 lc rgb "#808080"    t "{/Times-Italic t} = 10^{10}, inverse-time", \
    "../../data/scan/sigscan_t1e10_l.dat"  u 1:($6**2) w l lt 5 lw 4 lc rgb "#808080"    t "{/Times-Italic t} = 10^{10}, exact", \
    -1 notitle



unset multiplot
unset output
set terminal pop
reset
