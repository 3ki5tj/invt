#!/usr/bin/env gnuplot



# Normalized error vs the width of the Gaussian updating scheme



set encoding cp1250 # make the minus sign longer
set terminal push
# dl 3.0 make dashed line longer
set terminal postscript eps enhanced dl 3.0 size 8, 4 font "Times, 28"
set output "alphacmp.eps"
set multiplot


reset

wleft = 0.53
wright = 1 - wleft

set size wleft, 1
set origin 0, 0

set xtics 5e7 offset 0, 0.3
set format x "%.0t{/Symbol \264}10^{%T}"
set mxtics 5
set xlabel "Time, {/Times-Italic t}" offset 0, 0.5 font "Times, 36"

set logscale y
set format y "10^{%T}"
set ylabel "{/Symbol-Oblique a} ({/Times-Italic t})" font "Times, 36"

# `width` to reduce the text length
set key Left reverse width -6 spacing 1.4

set title "Nearest-neighbor updating scheme, {/Symbol-Oblique m}_1 = 0.24"

plot [:][5e-9:1e-4] \
    "../../data/scan/nb0.24_g_alpha.dat"  u 1:($2) w l lt 1 lw 4 lc rgb "#000000"    t "Global, inverse-time", \
    "../../data/scan/nb0.24_g_alpha.dat"  u 1:($4) w l lt 2 lw 4 lc rgb "#000000"    t "Global, exact", \
    "../../data/scan/nb0.24_l_alpha.dat"  u 1:($2) w l lt 4 lw 6 lc rgb "#808080"    t "Local, inverse-time", \
    "../../data/scan/nb0.24_l_alpha.dat"  u 1:($4) w l lt 5 lw 4 lc rgb "#808080"    t "Local, exact", \
    -1 notitle




set origin wleft, 0
set size wright, 1

unset ylabel

# `width` to reduce the text length
set key Left reverse width -5 spacing 1.4

set title "Gaussian updating scheme, {/Symbol-Oblique s} = 5"

set lmargin 4

plot [:][5e-9:1e-4] \
    "../../data/scan/sig5_g_alpha.dat"  u 1:($2) w l lt 1 lw 4 lc rgb "#000000"    t "Global, inverse-time", \
    "../../data/scan/sig5_g_alpha.dat"  u 1:($4) w l lt 2 lw 4 lc rgb "#000000"    t "Global, exact", \
    "../../data/scan/sig5_l_alpha.dat"  u 1:($2) w l lt 4 lw 6 lc rgb "#808080"    t "Local, inverse-time", \
    "../../data/scan/sig5_l_alpha.dat"  u 1:($4) w l lt 5 lw 4 lc rgb "#808080"    t "Local, exact", \
    -1 notitle



unset multiplot
unset output
set terminal pop
reset
