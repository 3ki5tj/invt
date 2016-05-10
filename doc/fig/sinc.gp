#!/usr/bin/env gnuplot



# optimized schedule and the bandpass schemes
# for Gaussian updating scheme with sigma = 5



set encoding cp1250 # make the minus sign longer
set terminal push
# dl 3 make dashed line longer
set terminal postscript eps enhanced dl 3 size 3.7, 3.7 font "Times, 24"
set output "sinc.eps"
set multiplot


dx = 0.01
dy = 0.05

normd(x, sig) = exp(-x*x/2/sig/sig)/sqrt(2*pi)/sig;

# bottom panel
set size 1, 1
set origin 0, 0

set tmargin 0.5
set lmargin 7

set xtics offset 0, 0.2
set mxtics 10
set xlabel "{/Times-Italic i}" offset 0, 0.5

set ytics 0.02 offset 0.2, 0
set mytics 2
set ylabel "{/Symbol-Oblique m}_{/Times-Italic i}" offset 2, -1

set key right bottom Left reverse invert spacing 1.2 width -8

plot [:100][-0.04:0.07] \
    0 lt 1 lw 0.1 lc rgb "#cccccc" notitle, \
    normd(x, 10) lt 1 lw 5 lc rgb "#aaaaaa" t "Gaussian, {/Symbol-Oblique s} = 10", \
    "../../data/sinc/sinc_nonpbc_win.dat" u 1:($2) w lp pt 4 ps 1.0 lt 1 lw 1.0 t "Bandpass, non{/*0.7 -}periodic, {/Times-Italic K} = 4", \
    "../../data/sinc/sinc_pbc_win.dat"    u 1:($2) w lp pt 5        lt 4 lw 1.0 t "Bandpass, periodic, {/Times-Italic K} = 2"


# inset
set origin 0.48, 0.51
set size 0.48, 0.57


set title "{/Times-Italic w_{ij}}" offset 0, -0.7 font "Times, 20"
#unset colorbox
set cbrange [-0.03:0.09]
set cbtics 0.05 offset -1, 0 font "Times, 14"
set mcbtics 5

set xtics 20 offset 0, 0.8 font "Times, 14"
set mxtics 2
#set format x ""
set xlabel "{/Times-Italic*0.7 j}" offset 0, 1.8

set ytics 20 offset 1, 0 font "Times, 14"
set mytics 2
#set format y ""
set ylabel "{/Times-Italic*0.7 i}" offset 3.2, 0

set tmargin 0
set bmargin 0
set lmargin 0
set rmargin 0

set pm3d map
splot [:] "../../data/sinc/sinc_nonpbc_winmat.dat" u 2:1:3 notitle

unset pm3d

unset multiplot
unset output
set terminal pop
reset
