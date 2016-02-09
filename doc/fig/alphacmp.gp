#!/usr/bin/env gnuplot



# Normalized error vs the width of the Gaussian updating scheme



set encoding cp1250 # make the minus sign longer
set terminal push
# dl 3.0 make dashed line longer
set terminal postscript eps enhanced dl 3.0 size 3.5, 5 font "Times, 24"
set output "alphacmp.eps"
set multiplot


reset

htop = 0.48
hbot = 1 - htop

# bottom panel
set size 1, hbot
set origin 0, 0

#set xtics 5e7 offset 0, 0.3
#set format x "%.0t{/Symbol \264}10^{%T}"
set xtics ("0" 0, \
           "" 1e7 1, \
           "" 2e7 1, \
           "" 3e7 1, \
           "" 4e7 1, \
           "5{/Symbol \264}10^7" 5e7, \
           "" 6e7 1, \
           "" 7e7 1, \
           "" 8e7 1, \
           "" 9e7 1, \
           "10^8" 1e8) offset 0, 0.2
set mxtics 5
set xlabel "Time, {/Times-Italic t}" offset 0, 0.5

set logscale y
set format y "10^{%T}"
set mytics 10
set ylabel "{/Symbol-Oblique a} ({/Times-Italic t})"

# `width` to reduce the text length
set key Left reverse width -6 spacing 1.2

set title "Gaussian updating scheme, {/Symbol-Oblique s} = 5" font "Times, 22"

plot [:1e8][8e-9:3e-4] \
    "../../data/scan/sig5_g_alpha.dat"  u 1:($2) w l lt 1 lw 4 lc rgb "#000000"    t "Global, inverse-time", \
    "../../data/scan/sig5_g_alpha.dat"  u 1:($4) w l lt 2 lw 4 lc rgb "#000000"    t "Global, exact", \
    "../../data/scan/sig5_l_alpha.dat"  u 1:($2) w l lt 4 lw 6 lc rgb "#808080"    t "Local, inverse-time", \
    "../../data/scan/sig5_l_alpha.dat"  u 1:($4) w l lt 5 lw 4 lc rgb "#808080"    t "Local, exact", \
    -1 notitle



# top panel
set origin 0, hbot
set size 1, htop

# the top panel share the same x-axis with
# the bottom one
unset xlabel
set format x ""



set title "Nearest-neighbor updating scheme, {/Symbol-Oblique m}_1 = 0.24"

plot [:1e8][8e-9:1e-4] \
    "../../data/scan/nb0.24_g_alpha.dat"  u 1:($2) w l lt 1 lw 4 lc rgb "#000000"    t "Global, inverse-time", \
    "../../data/scan/nb0.24_g_alpha.dat"  u 1:($4) w l lt 2 lw 4 lc rgb "#000000"    t "Global, exact", \
    "../../data/scan/nb0.24_l_alpha.dat"  u 1:($2) w l lt 4 lw 6 lc rgb "#808080"    t "Local, inverse-time", \
    "../../data/scan/nb0.24_l_alpha.dat"  u 1:($4) w l lt 5 lw 4 lc rgb "#808080"    t "Local, exact", \
    -1 notitle



unset multiplot
unset output
set terminal pop
reset
