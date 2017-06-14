#!/usr/bin/env gnuplot



# Normalized error vs the width of the Gaussian updating scheme



set encoding cp1250 # make the minus sign longer
set terminal push
# dl 3.0 make dashed line longer
set terminal postscript eps enhanced dl 3.0 size 3.5, 5 font "Times, 24"
set output "alphacmp.eps"
set multiplot


reset

htop = 0.49
hbot = 1 - htop

# bottom panel
set size 1, hbot
set origin 0, 0

set format x "%.0t{/Symbol \264}10^{%T}"
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
#set logscale x
#set mxtics 10
#set format x "10^{%T}"
set xrange [1e4:1e8]
set xlabel "Time, {/Times-Italic t}" offset 0, 0.5

#set logscale y
#set format y "10^{%T}"
#set format y "%.0t{/Symbol \264}10^{%T}"
set ytics ("0" 0, \
           "" 1e7 1, \
           "" 2e7 1, \
           "" 3e7 1, \
           "" 4e7 1, \
           "5{/Symbol \264}10^7" 5e7, \
           "" 6e7 1, \
           "" 7e7 1, \
           "" 8e7 1, \
           "" 9e7 1, \
           "10^8" 1e8, \
           "" 11e7 1, \
           "" 12e7 1 \
) offset 0.2, 0
set mytics 5
set yrange [0:1.2e8]
set ylabel "1 / {/Symbol-Oblique a} ({/Times-Italic t})"

set title "Gaussian updating scheme, {/Symbol-Oblique s} = 5" offset -3, -0.2

set key left top Left reverse spacing 1.1

plot [][] \
    "../../data/scan/sig5_g_alpha.dat"  u 1:(1/$2) w l lt 1 lw 6 lc rgb "#000000"    t "Global, optimal", \
    "../../data/scan/sig5_g_alpha.dat"  u 1:(1/$4) w l lt 2 lw 6 lc rgb "#000000"    t "Global, inv. {/Times-Italic t}", \
    "../../data/scan/sig5_l_alpha.dat"  u 1:(1/$2) w l lt 1 lw 2 lc rgb "#000000"    t "Local, optimal", \
    "../../data/scan/sig5_l_alpha.dat"  u 1:(1/$4) w l lt 2 lw 2 lc rgb "#000000"    t "Local, inv. {/Times-Italic t}", \
    -1 notitle



# top panel
set origin 0, hbot
set size 1, htop

# the top panel share the same x-axis with
# the bottom one
unset xlabel
#set format x ""


set title "Nearest-neighbor updating scheme, {/Symbol-Oblique m}_1 = 0.24"

plot [][] \
    "../../data/scan/nb0.24_g_alpha.dat"  u 1:(1/$2) w l lt 1 lw 6 lc rgb "#000000"    t "Global, optimal", \
    "../../data/scan/nb0.24_g_alpha.dat"  u 1:(1/$4) w l lt 2 lw 6 lc rgb "#000000"    t "Global, inv. {/Times-Italic t}", \
    "../../data/scan/nb0.24_l_alpha.dat"  u 1:(1/$2) w l lt 1 lw 2 lc rgb "#000000"    t "Local, optimal", \
    "../../data/scan/nb0.24_l_alpha.dat"  u 1:(1/$4) w l lt 2 lw 2 lc rgb "#000000"    t "Local, inv. {/Times-Italic t}", \
    -1 notitle



unset multiplot
unset output
set terminal pop
reset
