#!/usr/bin/env gnuplot



# Normalized error vs the width of the Gaussian updating scheme



set encoding cp1250 # make the minus sign longer
set terminal push
# dl 3.0 make dashed line longer
set terminal postscript eps enhanced dl 3.0 size 3.5, 5 font "Times, 24"
set output "dinvacmp.eps"
set multiplot


reset

htop = 0.45
hbot = 1 - htop

# bottom panel
set size 1, hbot
set origin 0, 0

set logscale x
set mxtics 10
set format x "10^{%T}"
set xrange [3e4:1e8]
set xlabel "Time, {/Times-Italic t}" offset 0, 0.0

tmin = 3e4
tmax = 9.9e7

set logscale y
set format y "10^{%T}"
set mytics 10
set yrange [:20]
set ylabel "{/Times-Italic d} [1/{/Symbol-Oblique a}({/Times-Italic t})] / {/Times-Italic d t}"



set title "Gaussian updating scheme, {/Symbol-Oblique s} = 5" offset -1.0, -0.2

set key left top Left reverse spacing 1.1

sig_c_g = `tail -n 1 ../../data/scan/sig5_g_alpha.dat | cut -f 7`
sig_c_l = `tail -n 1 ../../data/scan/sig5_l_alpha.dat | cut -f 7`

plot [][:2e1] \
    "../../data/scan/sig5_g_alpha.dat"  u 1:(($1 > tmin && $1 < tmax) ? $4 : 1/0) w l lt 1 lw 6 t "Global, optimal", \
    sig_c_g                                                                       w l lt 2 lw 6 t "Global, inv. {/Times-Italic t}", \
    "../../data/scan/sig5_l_alpha.dat"  u 1:(($1 > tmin && $1 < tmax) ? $4 : 1/0) w l lt 1 lw 2 t "Local, optimal", \
    sig_c_l                                                                       w l lt 2 lw 2 t "Local, inv. {/Times-Italic t}", \
    -1 notitle



# top panel
set origin 0, hbot
set size 1, htop

set bmargin 0

# the top panel share the same x-axis with
# the bottom one
unset xlabel
set format x ""

nb_c_g = `tail -n 1 ../../data/scan/nb0.24_g_alpha.dat | cut -f 7`
nb_c_l = `tail -n 1 ../../data/scan/nb0.24_l_alpha.dat | cut -f 7`

set title "Nearest-neighbor updating scheme, {/Symbol-Oblique m}_1 = 0.24"

plot [][4e-2:2e1] \
    "../../data/scan/nb0.24_g_alpha.dat"  u 1:(($1 > tmin && $1 < tmax) ? $4 : 1/0) w l lt 1 lw 6 t "Global, optimal", \
    nb_c_g                                                                          w l lt 2 lw 6 t "Global, inv. {/Times-Italic t}", \
    "../../data/scan/nb0.24_l_alpha.dat"  u 1:(($1 > tmin && $1 < tmax) ? $4 : 1/0) w l lt 1 lw 2 t "Local, optimal", \
    nb_c_l                                                                          w l lt 2 lw 2 t "Local, inv. {/Times-Italic t}", \
    -1 notitle



unset multiplot
unset output
set terminal pop
reset
