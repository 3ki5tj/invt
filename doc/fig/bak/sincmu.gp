#!/usr/bin/env gnuplot



# optimized schedule of the bandpass updating schemes



set encoding cp1250 # make the minus sign longer
set terminal push
# dl 4 make dashed line longer
set terminal postscript eps enhanced dl 4 size 3.5, 5.5 font "Times, 24"
set output "sincmu.eps"
set multiplot


reset

htop = 0.35
hbot = 1 - htop

dx = 0.02
dy = 0.05

set label "(a)" at screen dx, 1 - dy    font "Times, 32"
set label "(b)" at screen dx, hbot - dy font "Times, 32"

normd(x, sig) = exp(-x*x/2/sig/sig)/sqrt(2*pi)/sig;

# bottom panel
set size 1, hbot
set origin 0, 0

set tmargin 0
#set lmargin 8

set pm3d map

set xtics 20
set mxtics 2

set ytics 20
set mytics 2

splot [][] \
    "../../data/scan/sinc_nonpbc_winmat.dat" u 2:1:3 notitle

unset pm3d

reset

# top panel
set origin 0, hbot
set size 1, htop

set bmargin 0
set xtics offset 0, 0.2
set mxtics 10
set xlabel "{/Times-Italic i}" offset 0, 0.5

set ytics 0.05 offset 0.2, 0
set mytics 5
set ylabel "{/Symbol-Oblique m}_{/Times-Italic i}" offset 2, -1

set key invert spacing 1.2

plot [:100][-0.03:0.1] \
    0 lt 1 lw 0.1 lc rgb "#aaaaaa" notitle, \
    normd(x, 5) lt 1 lw 4 lc rgb "#cccccc" t "Gaussian, {/Symbol-Oblique s} = 5", \
    "../../data/scan/sinc_nonpbc_win.dat" u 1:($2) w lp pt 2 ps 0.8 lt 1 lw 1 lc rgb "#000080" t "Non{/*0.7 -}periodic, {/Times-Italic K} = 9", \
    "../../data/scan/sinc_pbc_win.dat"    u 1:($2) w lp pt 6        lt 1 lw 1 lc rgb "#800000" t "Periodic, {/Times-Italic K} = 4"


unset multiplot
unset output
set terminal pop
reset
