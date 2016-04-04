#!/usr/bin/env gnuplot



# optimized schedule versus the optimal inverse-time schedule
# for Gaussian updating scheme with sigma = 5



set encoding cp1250 # make the minus sign longer
set terminal push
# dl 4 make dashed line longer
set terminal postscript eps enhanced dl 4 size 3.5, 4.5 font "Times, 24"
set output "sinc.eps"
set multiplot


reset

htop = 0.5
hbot = 1 - htop

dx = 0.02
dy = 0.05

set label "(a)" at screen dx, 1 - dy    font "Times, 32"
set label "(b)" at screen dx, hbot - dy font "Times, 32"

normd(x, sig) = exp(-x*x/2/sig/sig)/sqrt(2*pi)/sig;

# bottom panel
set size 1, hbot
set origin 0, 0

set lmargin 7
set xlabel "{/Times-Italic i}"
set mxtics 10

set ytics 0.05
set mytics 5
set ylabel "{/Symbol-Oblique m}_{/Times-Italic i}" offset 3, 0

plot [:50][-0.03:0.1] \
    "../../data/scan/sinc_pbc_win.dat"    u 1:($2) w lp pt 6 lt 3 lw 2 t "Periodic, {/Times-Italic K} = 4", \
    "../../data/scan/sinc_nonpbc_win.dat" u 1:($2) w lp pt 2 lt 1 lw 2 t "Nonperiodic, {/Times-Italic K} = 8", \
    normd(x, 4.5) lt 2 lw 5 lc rgb "#808080" t "Gaussian, {/Symbol-Oblique s} = 4.5", \
    0 lt 1 lw 0.1 lc rgb "#aaaaaa" notitle

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

plot [][2e-9:5e-5] \
    1/(x + 20000) lw 2 t "1 / ({/Times-Italic t} + 2/{/Symbol-Oblique a}_0)", \
    "../../data/scan/sinc_pbc_alpha.dat"    u 1:($2) every 300 w p pt 6 ps 1.4 lw 2 t "Periodic, {/Times-Italic K} = 4", \
    "../../data/scan/sinc_nonpbc_alpha.dat" u 1:($2) every 300 w p pt 2 ps 1.4 lw 2 t "Nonperiodic, {/Times-Italic K} = 8", \
    -1 notitle




unset multiplot
unset output
set terminal pop
reset
