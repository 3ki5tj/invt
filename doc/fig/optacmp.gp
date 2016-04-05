#!/usr/bin/env gnuplot



# optimized schedule versus the optimal inverse-time schedule
# for Gaussian updating scheme with sigma = 5



set encoding cp1250 # make the minus sign longer
set terminal push
# dl 2.0 make dashed line longer
set terminal postscript eps enhanced dl 2.0 size 3.5, 4.5 font "Times, 24"
set output "optacmp.eps"
set multiplot


reset

htop = 0.45
hbot = 1 - htop

dx = 0.02
dy = 0.05

set label "(a)" at screen dx, 1 - dy    font "Times, 32"
set label "(b)" at screen dx, hbot - dy font "Times, 32"

# bottom panel
set size 1, hbot
set origin 0, 0

set logscale x
set mxtics 10
set format x "10^{%T}"
set xrange [1e4:1e8]
set xlabel "Simulation time, {/Times-Italic t}" offset 0, 0.0

set logscale y
set format y "10^{%T}"
set mytics 10
#set yrange [1e4:1e9]
set ylabel "{/Times-Italic d} {/Symbol-Oblique a}^{/*0.7 -1}({/Times-Italic t}) / {/Times-Italic d t}"

set key left Left reverse

sig_c_g = `tail -n 1 ../../data/scan/sig10_g_alpha.dat | cut -f 7`
sig_c_l = `tail -n 1 ../../data/scan/sig10_l_alpha.dat | cut -f 7`

plot [][:20] \
    "../../data/scan/sig10_g_alpha.dat"  u 1:4 w l lt 1 lw 6 t "Global, optimal", \
    sig_c_g                                    w l lt 2 lw 6 t "Global, inverse {/Times-Italic t}", \
    "../../data/scan/sig10_l_alpha.dat"  u 1:4 w l lt 1 lw 2 t "Local, optimal", \
    sig_c_l                                    w l lt 2 lw 2 t "Local, inverse {/Times-Italic t}", \
    -1 notitle


# top panel
set origin 0, hbot
set size 1, htop

set bmargin 0

# the top panel share the same x-axis with
# the bottom one
unset xlabel
set format x ""

set ylabel "{/Symbol-Oblique a} ({/Times-Italic t})"


set key left bottom Left reverse spacing 1.1

plot [][2e-9:2e-4] \
    "../../data/scan/sig10_g_alpha.dat"  u 1:($2) w l lt 1 lw 6 lc rgb "#000000"    t "Global, optimal", \
    "../../data/scan/sig10_g_alpha.dat"  u 1:($5) w l lt 2 lw 6 lc rgb "#000000"    t "Global, inverse {/Times-Italic t}", \
    "../../data/scan/sig10_l_alpha.dat"  u 1:($2) w l lt 1 lw 2 lc rgb "#000000"    t "Local, optimal", \
    "../../data/scan/sig10_l_alpha.dat"  u 1:($5) w l lt 2 lw 2 lc rgb "#000000"    t "Local, inverse {/Times-Italic t}", \
    -1 notitle




unset multiplot
unset output
set terminal pop
reset
