#!/usr/bin/env gnuplot



# optimized schedule versus the optimal inverse-time schedule
# for Gaussian updating scheme with sigma = 5



set encoding cp1250 # make the minus sign longer
set terminal push
# dl 2.0 make dashed line longer
set terminal postscript eps enhanced dl 2.0 size 3.5, 5.0 font "Times, 24"
set output "optacmp.eps"
set multiplot


reset

htop = 0.5
hbot = 1 - htop

dx = 0.01
dy = 0.05

set label "(a)" at screen dx, 1 - dy    font "Times, 32"
set label "(b)" at screen dx, hbot - dy font "Times, 32"

# bottom panel
set size 1, hbot
set origin 0, 0

set xtics 10
set mxtics 10
set xlabel "{/Times-Italic q}({/Times-Italic T}) - {/Times-Italic q}"

set logscale y
set format y "10^{%T}"
set ylabel "{/Times-Italic m}({/Times-Italic q}) = 1/({/Times-Italic T}{/Symbol-Oblique a})"

set key right top Left reverse width -3 spacing 1.1

plot [:30][1e-3:] \
    "../../data/opta/sig10_g_alpha.dat"  u (494.499-$3):(1/(1e8*$2))           w l lt 1 lw 6  t "Gaussian, global", \
    "../../data/opta/sig10_l_alpha.dat"  u (162.972-$3):(1/(1e8*$2))           w l lt 1 lw 2  t "Gaussian, local", \
    "../../data/opta/wl_g_alpha.dat"     u (8.51733-$3):(1/(1e8*$2))           w l lt 4 lw 2  t "Single-bin, global", \
    "../../data/opta/wl_l_alpha.dat"     u (8.51733-$3):(1/(1e8*$2)) every 500 w p pt 6 lw 2  t "Single-bin, local", \
    0 notitle




reset
# top panel
set origin 0, hbot
set size 1, htop

set logscale x
set mxtics 10
set format x "10^{%T}"
set xrange [1e4:1e8]
set xlabel "Simulation time, {/Times-Italic t}" offset 0, 0.0

set logscale y
set format y "10^{%T}"
set mytics 10
#set yrange [1e4:1e9]
set ylabel "{/Symbol-Oblique a} ({/Times-Italic t})"


set key left bottom Left reverse spacing 1.1

plot [][3e-9:1e-4] \
    "../../data/opta/sig10_g_alpha.dat"  u 1:($2) w l lt 1 lw 6 lc rgb "#000000"    t "Global, optimal", \
    "../../data/opta/sig10_g_alpha.dat"  u 1:($5) w l lt 2 lw 6 lc rgb "#000000"    t "Global, inverse {/Times-Italic t}", \
    "../../data/opta/sig10_l_alpha.dat"  u 1:($2) w l lt 1 lw 2 lc rgb "#000000"    t "Local, optimal", \
    "../../data/opta/sig10_l_alpha.dat"  u 1:($5) w l lt 2 lw 2 lc rgb "#000000"    t "Local, inverse {/Times-Italic t}", \
    -1 notitle




unset multiplot
unset output
set terminal pop
reset
